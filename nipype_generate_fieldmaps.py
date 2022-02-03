from __future__ import annotations

import argparse
import json
import math
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

__version__ = "0.2.4"


INPUT_FIELDS = [
    "se_epi_pe1_file",
    "se_epi_pe2_file",
    "se_epi_pe1_sidecar_file",
    "se_epi_pe2_sidecar_file",
]
OUTPUT_FIELDS = [
    "acq_params_file",
    "echo_spacing",
    "pe1_pedir",
    "pe2_pedir",
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

PE_XYZ = {"i": "x", "j": "y", "k": "z", "i-": "-x", "j-": "-y", "k-": "-z"}


def create_generate_fieldmaps_wf(name: str = "generate_fieldmaps_wf") -> Workflow:

    wf = Workflow(name=name)
    wf.config["execution"]["remove_unnecessary_outputs"] = "false"

    inputnode = Node(IdentityInterface(fields=INPUT_FIELDS), name="inputnode")

    # extract the effective echo_spacing
    effective_echo_spacing = Node(
        Function(
            input_names=["sidecar_file"],
            output_names=["echo_spacing"],
            function=_extract_effective_echo_spacing_fi,
        ),
        name="effective_echo_spacing",
    )
    wf.connect(
        inputnode,
        "se_epi_pe1_sidecar_file",
        effective_echo_spacing,
        "sidecar_file",
    )

    # extract the epi_reg-compatible phase-encoding directions
    pe1_pedir = Node(
        Function(
            input_names=["sidecar_file"],
            output_names=["pedir"],
            function=_extract_pedir_fi,
        ),
        name="pe1_pedir",
    )
    wf.connect(inputnode, "se_epi_pe1_sidecar_file", pe1_pedir, "sidecar_file")

    pe2_pedir = Node(
        Function(
            input_names=["sidecar_file"],
            output_names=["pedir"],
            function=_extract_pedir_fi,
        ),
        name="pe2_pedir",
    )
    wf.connect(inputnode, "se_epi_pe2_sidecar_file", pe2_pedir, "sidecar_file")

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
                "pe1_sidecar_file",
                "pe2_sidecar_file",
                "out_file",
            ],
            output_names=["out_file"],
            function=_create_acq_param_file_fi,
        ),
        name="acq_params",
    )
    wf.connect(merge_se_epi_files, "merged_file", acq_params, "merged_se_epi_file")
    wf.connect(inputnode, "se_epi_pe1_sidecar_file", acq_params, "pe1_sidecar_file")
    wf.connect(inputnode, "se_epi_pe2_sidecar_file", acq_params, "pe2_sidecar_file")

    # estimate the fieldmaps via FSL's TOPUP
    topup = Node(
        fsl.TOPUP(out_field="fmap_hz.nii.gz", out_corrected="corrected.nii.gz"),
        name="topup",
    )
    wf.connect(merge_se_epi_files, "merged_file", topup, "in_file")
    wf.connect(acq_params, "out_file", topup, "encoding_file")

    # convert the estimated field to rad/s
    two_pi = 2 * math.pi
    fmap_rads = Node(
        fsl.ImageMaths(op_string=f"-mul {two_pi}", out_file="fmap_rads.nii.gz"),
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
    wf.connect(effective_echo_spacing, "echo_spacing", outputnode, "echo_spacing")
    wf.connect(pe1_pedir, "pedir", outputnode, "pe1_pedir")
    wf.connect(pe2_pedir, "pedir", outputnode, "pe2_pedir")
    wf.connect(topup, "out_corrected", outputnode, "corrected_se_epi_file")
    wf.connect(topup, "out_field", outputnode, "fmap_hz_file")
    wf.connect(fmap_rads, "out_file", outputnode, "fmap_rads_file")
    wf.connect(fmap_mag, "out_file", outputnode, "fmap_mag_file")
    wf.connect(fmap_mag_brain, "out_file", outputnode, "fmap_mag_brain_file")
    wf.connect(fmap_mag_brain, "mask_file", outputnode, "fmap_mag_brain_mask_file")

    return wf


# INTERFACES (fi = [F]unction [I]nterface)


def _extract_effective_echo_spacing_fi(sidecar_file):
    import json
    from pathlib import Path
    from nipype_generate_fieldmaps import get_effective_echo_spacing

    sidecar = json.loads(Path(sidecar_file).read_text())
    return get_effective_echo_spacing(sidecar)


def _extract_pedir_fi(sidecar_file):
    import json
    from pathlib import Path
    from nipype_generate_fieldmaps import get_phase_encoding_xyz

    sidecar = json.loads(Path(sidecar_file).read_text())
    return get_phase_encoding_xyz(sidecar)


def _create_acq_param_file_fi(
    merged_se_epi_file,
    pe1_sidecar_file,
    pe2_sidecar_file,
    out_file=None,
):
    from nipype_generate_fieldmaps import create_acq_param_file

    return create_acq_param_file(
        merged_se_epi_file,
        pe1_sidecar_file,
        pe2_sidecar_file,
        out_file,
    )


# UTILITIES


def create_acq_param_file(
    merged_se_epi_file: str | Path,
    pe1_sidecar_file: str | Path,
    pe2_sidecar_file: str | Path,
    out_file: str | Path | None = None,
) -> Path:
    # load JSON sidecars
    pe1_sidecar = json.loads(Path(pe1_sidecar_file).read_text())
    pe2_sidecar = json.loads(Path(pe2_sidecar_file).read_text())

    # total readout times
    trt_pe1 = get_total_readout_time(pe1_sidecar)
    trt_pe2 = get_total_readout_time(pe2_sidecar)

    # phase encoding unit vectors
    pe1_vec = get_phase_encoding_vec(pe1_sidecar)
    pe2_vec = get_phase_encoding_vec(pe2_sidecar)

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
    if "EffectiveEcho_Spacing" in sidecar:
        # can we extract it?
        effective_echo_spacing: float = sidecar["EffectiveEcho_Spacing"]
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
    return PE_XYZ[pe]


def cli() -> int:
    parser = create_parser()
    args = parser.parse_args()

    if hasattr(args, "handler"):
        return args.handler(args)

    parser.print_help()
    return 1


def create_parser(
    parser: argparse.ArgumentParser | None = None,
) -> argparse.ArgumentParser:
    description = (
        "Generate fieldmaps from EPI acquisitions with differing "
        "phase-encoding directions"
    )
    _parser = parser or argparse.ArgumentParser(description=description)
    _parser.add_argument("-v", "--version", action="version", version=__version__)
    _parser.add_argument(
        "se_epi_pe1",
        type=Path,
        help="The spin-echo EPI file acquired in the 'first' phase-encoding direction",
    )
    _parser.add_argument(
        "se_epi_pe2",
        type=Path,
        help="The spin-echo EPI file acquired in the 'second' phase-encoding direction",
    )
    _parser.add_argument(
        "se_epi_pe1_sidecar",
        type=Path,
        help="The JSON sidecar for the first spin-echo EPI file",
    )
    _parser.add_argument(
        "se_epi_pe2_sidecar",
        type=Path,
        help="The JSON sidecar for the second spin-echo EPI file",
    )
    _parser.add_argument(
        "out_dir",
        type=Path,
        help="The directory into which outputs are written",
    )

    _parser.set_defaults(handler=handler)

    return _parser


def handler(args: argparse.Namespace) -> int:
    se_epi_pe1: Path = args.se_epi_pe1
    se_epi_pe2: Path = args.se_epi_pe2
    se_epi_pe1_sidecar: Path = args.se_epi_pe1_sidecar
    se_epi_pe2_sidecar: Path = args.se_epi_pe2_sidecar
    out_dir: Path = args.out_dir

    # create the workflow
    wf = create_generate_fieldmaps_wf()

    # wire-up the workflow
    wf.base_dir = out_dir.expanduser().resolve()
    wf.inputs.inputnode.se_epi_pe1_file = se_epi_pe1.expanduser().resolve()
    wf.inputs.inputnode.se_epi_pe2_file = se_epi_pe2.expanduser().resolve()
    wf.inputs.inputnode.se_epi_pe1_sidecar_file = (
        se_epi_pe1_sidecar.expanduser().resolve()
    )
    wf.inputs.inputnode.se_epi_pe2_sidecar_file = (
        se_epi_pe2_sidecar.expanduser().resolve()
    )

    # run it!
    wf.run()

    return 0
