# Bibliography used in Paper A

Primary sources and their role (checked 3 October 2026):

- Zienkiewicz and Zhu (1992), [SPR, Part 1](https://doi.org/10.1002/nme.1620330702): classical recovery context. Our nodal cubic least-squares gradient control is ZZ-type, not the original superconvergent-point SPR construction.
- Naga and Zhang (2005), [higher-order polynomial-preserving recovery in 2D/3D](https://doi.org/10.3934/dcdsb.2005.5.769): polynomial preservation and recovery context, not an unconditional guarantee on our distorted meshes.
- Zhang (2004), [bilinear elements on irregular quadrilaterals](https://doi.org/10.4208/ijnam.OA-2004-1101): recovery under mesh assumptions.
- Ainsworth and Oden (2000), [A Posteriori Error Estimation in Finite Element Analysis](https://doi.org/10.1002/9781118032824), and Verfürth (2013), [A Posteriori Error Estimation Techniques for Finite Element Methods](https://doi.org/10.1093/acprof:oso/9780199679423.001.0001): residual estimation and reliability context. Our unit-constant control omits boundary-data oscillation and is not a certified implementation of those theories.
- Huang and Russell (2011), [Adaptive Moving Mesh Methods](https://doi.org/10.1007/978-1-4419-7916-2): equidistribution and alignment background.
- Dobrev et al. (2019), [high-order TMOP](https://doi.org/10.1137/18M1167206): target-matrix framework.
- Dobrev et al. (2020), [simulation-driven high-order mesh optimization](https://arxiv.org/abs/2001.11536): establishes why TMOP must not be described as necessarily solution-independent.
- Satheesh et al. (2024), [finite element distortion potentials](https://doi.org/10.1016/j.cma.2024.117444): regularization with localized requirements; cited, not implemented as an additional baseline here.
- [MFEM tmop.cpp](https://github.com/mfem/mfem/blob/master/fem/tmop.cpp), metric 002: the explicit shape metric 0.5 ||T||²/det(T)-1 independently implemented here. This source link is mutable; the paper and code state the actual formula, so reproduction does not depend on fetching it.

The manuscript's `references.bib` contains the corresponding entries and local references to Paper 1 and Paper B. Bibliographic context does not transfer convergence or reliability theorems to the current algorithm without their hypotheses.
