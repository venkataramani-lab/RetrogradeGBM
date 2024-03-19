# Atlas Registration: Defining distances of connected neurons to tumor

## Summary

To register and analyze brain sections, we used the [QUINT workflow](https://quint-workflow.readthedocs.io/en/latest/QUINTintro.html) consisting of three steps. First, the sections were registered to Allen Mouse Brain Common Coordinate Framework (CCF). Sections were then preprocessed and segmented for quantification.

## Workflow

1. Registration of brain sections to Allen Mouse Brain CCF using QUINT workflow
2. Segmentation of tumor using Fiji ImageJ
3. Segmentation of GFP positive cell centroids using Ilastik and Fiji ImageJ
4. Quantification using Nutil (part of QUINT workflow)
5. Using JSON coordinate files from Nutil: Matlab script "AtlasRegistration_DistanceToTumor.m" calculates distance of GFP positive cells to tumor

## Contact

For any questions, please feel free to reach out to us via email.