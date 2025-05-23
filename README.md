# Optimizing QAOA circuit transpilation with parity twine and SWAP network encodings
Mapping quantum approximate optimization algorithm (QAOA) circuits with non-trivial connectivity in fixed-layout quantum platforms such as superconducting-based quantum processing units (QPUs) requires a process of transpilation to match the quantum circuit on the given layout. This step is critical for reducing error rates when running on noisy QPUs. Two methodologies that improve the resource required to do such transpilation are the SWAP network and parity twine chains (PTC). These approaches reduce the two-qubit gate count and depth needed to represent fully connected circuits. In this work, a simulated annealing-based method is introduced that reduces the PTC and SWAP network encoding requirements in QAOA circuits with non-fully connected two-qubit gates. This method is benchmarked against various transpilers and demonstrates that, beyond specific connectivity thresholds, it achieves significant reductions in both two-qubit gate count and circuit depth, surpassing the performance of Qiskit’s transpiler at its highest optimization level. For example, for a 120-qubit QAOA instance with 25\% connectivity, our method achieves an 85\% reduction in depth and a 28\% reduction in two-qubit gates. Finally, the practical impact of PTC encoding is validated by benchmarking QAOA on the \texttt{ibm\_fez} device, showing improved performance up to 20 qubits, compared to a 15-qubit limit when using SWAP networks.


![Alt text](paper-transpilation.png)



