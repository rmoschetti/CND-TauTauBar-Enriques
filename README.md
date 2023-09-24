# CND-TauTauBar-Enriques
This project contains the computational data of the paper *"The non-degeneracy invariant of Brandhorst and Shimada families of Enriques surfaces"*, by Riccardo Moschetti, Franco Rota, and Luca Schaffler.
The data was obtained by applying the code *CNDFinder* [MRS22a, MRS22b] to Brandhorst and Shimada's data [SB20] about smooth rational curves on realizable $(\tau,\overline{\tau})$-generic Enriques surfaces, see [BS22].
As explained more in detail in the paper, for each realizable $(\tau,\overline{\tau})$-generic Enriques surfaces $Y$ we provide a non-degenerate isotropic sequence. The length of this sequence gives a lower bound for the non-degeneracy invariant of $Y$. In 141 of the 155 cases for $Y$ we obtain a sequence of length 10, which shows that the non-degeneracy invariant of $Y$ is 10. Isotropic vectors and smooth rational curves are expressed as vectors whose coordinates are with respect to a fixed $E_{10}$-basis for $\mathrm{Num}(Y)$, as done in [BS22].

**Summary of the files contained in this folder**
 - Cnd_Data.txt - a list with the computational data we obtained. The data are formatted as JSON text object, and can be imported in SageMath.
 - CheckData.py - a SageMath code which loads "Cnd_Data.txt" and automatically performs all the necessary checks to ensure that the sequences in "Cnd_Data.txt" are indeed non-degenerate isotropic sequences. This gives a computational proof of Theorem 4.2 (2) and (3).
 - Parser.py - we do not include the original data from Brandhorst and Shimada, which can be found in [SB20]. However, for the user convenience, this parser takes as input the original file [SB20] (which are made in GAP), and create the JSON files, extracting only the data of smooth rational curves and automorphisms for each case of realizable $(\tau,\overline{\tau})$-generic Enriques surface.

**Data structure of Cnd_Data.txt**
Cnd_Data.txt contains a list of objects, one for each realizable $(\tau,\overline{\tau})$-generic Enriques surface. For each object we have the following keys:

 - "**no"**: the number identifying the $(\tau,\overline{\tau})$-generic Enriques surface as listed in [Table 1, BS22].
 - "**cnd**": the lower bound for the non-degeneracy invariant that we computed.
 - "**type**": a string describing the type of elliptic configurations appearing in the non-degenerate isotropic sequence.
 - "**tex**": a string containing the isotropic sequence written in $\LaTeX$.
 - "**IsotropicSequence**": a dictionary describing the isotropic sequence for computational purposes (this is the "core" of these data). The dictionary contains a key "EllipticConfiguration" which is a list of lenght "cnd" of elliptic configurations, and a key "Coefficient", which is 1 or 1/2. Each isotropic vector in the sequence is described as in the following examples:

> {'Coefficient': 1, 'EllipticConfiguration': [[1, 'R0'], [1, 'R29']]}

The above gives the half-fiber with numerical equivalence class $(1\cdot R_{0} + 1 \cdot R_{29})$. The explicit values of R0 and R29 are taken from "RatsUsed".

> {'Coefficient': 1/2,    'EllipticConfiguration': [[1, 'R0'], [1, 'R8'], [1, 'R5', 'H2'], [1, 'R1', 'H4'], [2, 'R5']]}

The above gives the half-fiber with numerical equivalence class $\frac{1}{2}(1\cdot R_{0} + 1 \cdot R_{8} + 1 \cdot R_{5} \cdot H_{2} + 1 \cdot R_{1} \cdot H_{4} + 2 \cdot R_{5})$. The explicit values of R0, R8, R5, R1 are taken from "RatsUsed", while H2, H4 taken from "AutUsed".

> {'Coefficient': 1/2, 'EllipticConfiguration': [[1, 'R2', 'H6'], [1, 'R3', 'H6', 'H0']]}

The above gives the half-fiber with numerical equivalence class $\frac{1}{2}(1\cdot R_{2} \cdot H_{6} + 1 \cdot R_{3} \cdot H_{6} \cdot H_{0})$. The explicit values of R2, R3 are taken from "RatsUsed" while H0, H6 are taken from "AutUsed".

 - "**RatsUsed**": an object which describes a subset of the list "Ratstemp" in [SB20], containing only the coordinates of the smooth rational curves which were used in the description of the corresponding isotropic sequence. The curves are labeled RX, where X is the position in the list provided in [SB20], starting from zero.
 - "**AutUsed**": an object which describes a subset of the list "Autrec" in [SB20], containing only the $10 \times 10$ matrices of automorphisms which were used in the description of the corresponding isotropic sequence. *Note that this is empty for most of the cases.* The automorphisms are labeled HX, where X is the position in the list provided in [SB20], starting from zero.
 - "**ProducedIsotropicSequence**": this is the list of coordinates of the elliptic configurations listed in "IsotropicSequence". More precisely, if we have the numerical equivalence class of an half-fiber given by $\frac{1}{2}(1\cdot R_{2} \cdot H_{6} + 1 \cdot R_{3} \cdot H_{6} \cdot H_{0})$, here we expand such an expression and store the resulting length 10 vector.

To use the code, just open in SageMath a new project in the same folder as the file Cnd_Data.txt. Copy-paste the content of CheckData.py in the project, and run it.

[BS22] Brandhorst, Simon and Shimada, Ichiro; Automorphism groups of certain Enriques surfaces, Found. Comput. Math. 22, 2022.

[MRS22a] Riccardo Moschetti, Franco Rota, and Luca Schaffler. A computational view on the non-degeneracy invariant for Enriques surfaces. Experimental Mathematics. Published online: 29 August 2022.

[MRS22b] Riccardo Moschetti, Franco Rota, and Luca Schaffler. SageMath code CndFinder. Available at: https://github.com/rmoschetti/CNDFinder

[SB20] Shimada, Ichiro and Brandhorst, Simon; Automorphism groups of certain Enriques surfaces - computational data, https://zenodo.org/record/4327019