from __future__ import annotations

import json
from itertools import chain
from pathlib import Path
from typing import Any

import nibabel as nib
from nipype import Function
from nipype import IdentityInterface
from nipype import Merge
from nipype import Node
from nipype import Workflow
from nipype.interfaces import fsl

__version__ = "0.1.0"


INPUT_FIELDS = [
    "se_epi_pe1_file",
    "se_epi_pe2_file",
    "se_epi_sidecar_pe1_file",
    "se_epi_sidecar_pe2_file",
]
OUTPUT_FIELDS = [
    "acq_params_file",
    "echospacing",
    "pedir",
    "corrected_se_epi_file",
    "fmap_hz_file",
    "fmap_rads_file",
    "fmap_mag_file",
    "fmap_mag_brain_file",
    "fmap_mag_brain_mask_file",
]


PE_UVECTORS: dict[str, tuple[int, int, int]] = {  # unit vectors
    "i": (1, 0, 0),
    "j": (0, 1, 0),
    "k": (0, 0, 1),
    "i-": (-1, 0, 0),
    "j-": (0, -1, 0),
    "k-": (0, 0, -1),
}
EPI_REG_PEDIRS = {"i": "x", "j": "y", "k": "z", "i-": "-x", "j-": "-y", "k-": "-z"}


def create_prepare_fieldmaps_wf(name: str = "prepare_fieldmaps_wf") -> Workflow:

    wf = Workflow(name=name)
    wf.config["execution"]["remove_unnecessary_outputs"] = "false"

    inputnode = Node(IdentityInterface(fields=INPUT_FIELDS), name="inputnode")

    # extract the effective echospacing
    effective_echospacing = Node(
        Function(
            input_names=["sidecar_file"],
            output_names=["effective_echospacing"],
            function=_extract_effective_echo_spacing_fi,
        ),
        name="effective_echo_spacing",
    )
    wf.connect(
        inputnode,
        "se_epi_sidecar_pe1_file",
        effective_echospacing,
        "sidecar_file",
    )

    # extract the epi_reg-compatible phase-encoding direction
    pedir = Node(
        Function(
            input_names=["sidecar_file"],
            output_names=["pedir"],
            function=_extract_pedir_fi,
        ),
        name="pedir",
    )
    wf.connect(inputnode, "se_epi_sidecar_pe1_file", pedir, "sidecar_file")

    # pre-concatenation (need images in a list)
    listify_se_epi_files = Node(Merge(numinputs=2), name="listify_se_epi_files")
    wf.connect(inputnode, "se_epi_pe1_file", listify_se_epi_files, "in1")
    wf.connect(inputnode, "se_epi_pe2_file", listify_se_epi_files, "in2")

    # merge the acquisitions (volumes) into a single nii image
    merge_se_epi_files = Node(
        fsl.Merge(dimension="t", merged_file="merged.nii.gz"),
        name="merge_se_epi_files",
    )
    wf.connect(listify_se_epi_files, "out", merge_se_epi_files, "in_files")

    # create the acquisition parameter file (--datain)
    acq_params = Node(
        Function(
            input_names=[
                "merged_se_epi_file",
                "sidecar_pe1_file",
                "sidecar_pe2_file",
                "out_file",
            ],
            output_names=["out_file"],
            function=_create_acq_param_file_fi,
        ),
        name="acq_params",
    )
    wf.connect(merge_se_epi_files, "merged_file", acq_params, "merged_se_epi_file")
    wf.connect(inputnode, "se_epi_sidecar_pe1_file", acq_params, "sidecar_pe1_file")
    wf.connect(inputnode, "se_epi_sidecar_pe2_file", acq_params, "sidecar_pe2_file")

    # estimate the fieldmaps via FSL's TOPUP
    topup = Node(
        fsl.TOPUP(out_field="fmap_hz.nii.gz", out_corrected="corrected.nii.gz"),
        name="topup",
    )
    wf.connect(merge_se_epi_files, "merged_file", topup, "in_file")
    wf.connect(acq_params, "out_file", topup, "encoding_file")

    # convert the estimate field to rad/s
    fmap_rads = Node(
        fsl.ImageMaths(op_string="-mul 6.28", out_file="fmap_rads.nii.gz"),
        name="fmap_rads",
    )
    wf.connect(topup, "out_field", fmap_rads, "in_file")

    # compute a magnitude image from the corrected Spin Echo EPI volumes
    fmap_mag = Node(
        fsl.ImageMaths(op_string="-Tmean", out_file="fmap_mag.nii.gz"),
        name="fmap_mag",
    )
    wf.connect(topup, "out_corrected", fmap_mag, "in_file")

    # extract the mean brain + mask from the magnitude image
    fmap_mag_brain = Node(
        fsl.BET(frac=0.5, out_file="fmap_mag_brain.nii.gz", mask=True),
        name="fmap_mag_brain",
    )
    wf.connect(fmap_mag, "out_file", fmap_mag_brain, "in_file")

    # To the outside world!
    outputnode = Node(IdentityInterface(fields=OUTPUT_FIELDS), name="outputnode")
    wf.connect(acq_params, "out_file", outputnode, "acq_params_file")
    wf.connect(effective_echospacing, "echospacing", outputnode, "echospacing")
    wf.connect(pedir, "pedir", outputnode, "pedir")
    wf.connect(topup, "out_corrected", outputnode, "corrected_se_epi_file")
    wf.connect(topup, "out_field", outputnode, "fmap_hz_file")
    wf.connect(fmap_rads, "out_file", outputnode, "fmap_rads_file")
    wf.connect(fmap_mag, "out_file", outputnode, "fmap_mag_file")
    wf.connect(fmap_mag_brain, "out_file", outputnode, "fmap_mag_brain_file")
    wf.connect(fmap_mag_brain, "mask_file", outputnode, "fmap_mag_brain_mask_file")

    return wf


# INTERFACES (fi = [F]unction [I]nterface)


def _extract_effective_echo_spacing_fi(sidecar_file):
    from nipype_generate_fieldmaps import extract_effective_echo_spacing

    return extract_effective_echo_spacing(sidecar_file)


def _extract_pedir_fi(sidecar_file):
    from nipype_generate_fieldmaps import extract_pedir

    return extract_pedir(sidecar_file)


def _create_acq_param_file_fi(
    merged_se_epi_file,
    sidecar_pe1_file,
    sidecar_pe2_file,
    out_file=None,
):
    from nipype_generate_fieldmaps import create_acq_param_file

    return create_acq_param_file(
        merged_se_epi_file,
        sidecar_pe1_file,
        sidecar_pe2_file,
        out_file,
    )


# UTILITIES


def extract_effective_echo_spacing(sidecar_file: str | Path) -> float:
    sidecar = json.loads(Path(sidecar_file).read_text())
    return get_effective_echo_spacing(sidecar)


def extract_pedir(sidecar_file: str | Path) -> str:
    sidecar = json.loads(Path(sidecar_file).read_text())
    return get_phase_encoding_xyz(sidecar)


def create_acq_param_file(
    merged_se_epi_file: str | Path,
    sidecar_pe1_file: str | Path,
    sidecar_pe2_file: str | Path,
    out_file: str | Path | None = None,
) -> Path:
    # load JSON sidecars
    sidecar_pe1 = json.loads(Path(sidecar_pe1_file).read_text())
    sidecar_pe2 = json.loads(Path(sidecar_pe2_file).read_text())

    # total readout times
    trt_pe1 = get_total_readout_time(sidecar_pe1)
    trt_pe2 = get_total_readout_time(sidecar_pe2)

    # phase encoding unit vectors
    pe1_vec = get_phase_encoding_vec(sidecar_pe1)
    pe2_vec = get_phase_encoding_vec(sidecar_pe2)

    # extract the number of volumes in the merged fieldmaps nii image
    img: nib.Nifti1Image = nib.load(str(merged_se_epi_file))
    n_total_vols: int = img.header["dim"][4]  # type: ignore

    # format the lines that we'll write to the acq param file
    line_pe1 = " ".join(map(str, chain(pe1_vec, [trt_pe1])))
    lines_pe1 = [line_pe1] * (n_total_vols // 2)
    line_pe2 = " ".join(map(str, chain(pe2_vec, [trt_pe2])))
    lines_pe2 = [line_pe2] * (n_total_vols // 2)

    # create the acq param file
    content = "\n".join(chain(lines_pe1, lines_pe2)) + "\n"
    acq_param_file = Path(out_file) if out_file else Path.cwd() / "acq_params.txt"
    acq_param_file.write_text(content)

    return acq_param_file


def get_total_readout_time(sidecar: dict[str, Any]) -> float:
    # extract or derive the total readout time, see:
    # - https://bids-specification.readthedocs.io/en/v1.6.0/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#in-plane-spatial-encoding # noqa: E501
    # - https://lcni.uoregon.edu/kb-articles/kb-0003
    if "TotalReadoutTime" in sidecar:
        # can we extract it?
        total_readout_time: float = sidecar["TotalReadoutTime"]
    elif "ReconMatrixPE" in sidecar:
        # can we compute it?
        recon_matrix_pe: float = sidecar["ReconMatrixPE"]
        effective_echo_spacing = get_effective_echo_spacing(sidecar)
        total_readout_time = (recon_matrix_pe - 1) * effective_echo_spacing
    else:
        msg = "Could not extract or derive Total Readout Time from fieldmap sidecar"
        raise ValueError(msg)

    return total_readout_time


def get_effective_echo_spacing(sidecar: dict[str, Any]) -> float:
    if "EffectiveEchoSpacing" in sidecar:
        # can we extract it?
        effective_echo_spacing: float = sidecar["EffectiveEchoSpacing"]
    elif "BandwidthPerPixelPhaseEncode" in sidecar and "ReconMatrixPE" in sidecar:
        # can we compute it?
        bpppe: float = sidecar["BandwidthPerPixelPhaseEncode"]
        recon_matrix_pe: float = sidecar["ReconMatrixPE"]
        effective_echo_spacing = 1 / (bpppe * recon_matrix_pe)
    else:
        msg = "Could not extract or derive Effective Echo Spacing from fieldmap sidecar"
        raise ValueError(msg)

    return effective_echo_spacing


def get_phase_encoding_vec(sidecar: dict[str, Any]) -> tuple[int, int, int]:
    pe: str = sidecar["PhaseEncodingDirection"]
    return PE_UVECTORS[pe]


def get_phase_encoding_xyz(sidecar: dict[str, Any]) -> str:
    pe: str = sidecar["PhaseEncodingDirection"]
    return EPI_REG_PEDIRS[pe]
