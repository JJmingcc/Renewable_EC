# Optimal Workload Allocation for Distributed Edge Clouds With Renewable Energy and Battery Storage

This repository contains the official implementation of the proposed model as well as its compared benchmark. 

Pointers: [arxiv](https://arxiv.org/abs/2310.00742) [ICNC 2024](http://www.conf-icnc.org/2024/)


## Get started
### Prerequisites
- CVX: Matlab Software for Disciplined Convex Programming: [CXV](http://cvxr.com/cvx/)
- Gurobi solver (Free license for academic students): [Gurobi solver](https://www.gurobi.com/academia/academic-program-and-licenses/)
  
### File
- Model_1.m - Model_3.m: Models \textbf{M1} - \textbf{M3} encompass all the benchmarks we will use for comparison. These models correspond to those detailed in the Appendix of our report.
- Model_0.m: proposed model with battery and ``sell-back-to-grid" option
- Model_1.m: This model lacks the capability for the ``sell-back-to-grid" option, and ECs do not have %come equipped with batteries.
- Model_2.m: This model exclusively focuses on using batteries to store energy and does not allow the operator to sell surplus energy back to the grid.
- Model_3.m: This model is designed to enable the selling of excess energy back to the grid, while ECs do not feature battery installations.
- ICNC_results. mat: all saved data.
- results.zip: All saved figures used in paper
- init.m: initialization for parameter setting


## Citation and Acknowledgements
**Bibtex.**
If you find our code useful for your research, please cite the [arxiv](https://arxiv.org/abs/2310.00742):
```bibtex
@article{nguyen2023optimal,
  title={Optimal Workload Allocation for Distributed Edge Clouds With Renewable Energy and Battery Storage},
  author={Nguyen, Duong Thuy Anh and Cheng, Jiaming and Wang, Lele and Nguyen, Duong Tung},
  journal={arXiv preprint arXiv:2310.00742},
  year={2023}
}
```
**Acknowledgments and Disclosure of Funding.**
This work was funded, in part, by Natural Sciences and Engineering Research Council of Canada
(NSERC), Canada and National Science Foundation (NSF), USA.

## Contact
Please submit a GitHub issue or contact [jiaming@ece.ubc.ca](mailto:jiaming@ece.ubc.ca) if you have any questions or find any bugs.
