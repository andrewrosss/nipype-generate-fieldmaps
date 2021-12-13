# nipype-generate-fieldmaps

Nipype workflow to generate fieldmaps from EPI acquisitions with differing phase-encoding directions

## Installation

...

## Usage

...

## Prerequisites

- have two fieldmap acquisitions with differing phase encodings
- num vols in PE dir 1 == num vols in PE dir 2
- each fieldmap file has a sidecar with either:
  - TotalReadoutTime & PhaseEncodingDirection
  - EffectiveEchoSpacing & ReconMatrixPE & PhaseEncodingDirection
- need `graphviz` installed (potentially)
