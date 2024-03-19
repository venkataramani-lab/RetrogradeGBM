# Cluster Analysis of input and starter cells over time

## Summary

In this methodology, we employed a systematic analysis of input and starter cells and their spatial relationships. First, GFP and Nestin signals were segmented using ilastik. Starter cells were calculated by overlapping the segmented GFP and Nestin signals. Input cells were identified by subtracting starter cells from the GFP signal. We extracted the coordinates of the input and starter cells. This Matlab script creates clusters using the Density-Based Spatial Clustering of Applications with Noise (DBScan) algorithm. Parameters such as MinPts and epsilon can be adjusted according to the characteristics of individual samples. The MATLAB script calculates distances between starter cells, input cells, and the resulting clusters. Moreover, it transforms the cluster boundaries into Regions of Interest (ROI) represented by Convex Hulls for enhanced delineation.

## Contact

For any questions, please feel free to reach out to us via email.