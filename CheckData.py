from sage.all import *

def GetExplicitIsotropicSequence(d):
    # This function computes explicitly the coordinates of the isotropic vectors in d["IsotropicSequence"]
    # The isotropic vectors are obtained by combining smooth rational curves in d["RatsUsed"] and automorphisms in d["AutUsed"], together with some coefficients.
    # Explicitly, if an element in d["IsotropicSequence"] is 
    # {"Coefficient": "1/2", "EllipticConfiguration": [[1, "R2"], [1, "R14", "H0"]]}
    # then we first run the function GetCoefficientWithAutomorphisms to obtain the coordinates of R2
    # and the coordinates of R14 * H0, extracting the data from d["RatsUsed"] and d["AutUsed"].
    # After this, we compute 1/2(1*R2 + 1 *R14*H0) returning the corresponding length 10 integral vector.
    
    ProducedIsotropicSequence=[]
    for ec in d["IsotropicSequence"]:
        ExplicitVector=matrix(ZZ,[0]*10)
        for x in ec["EllipticConfiguration"]:
            ExplicitVector = ExplicitVector + Rational(x[0])*matrix(QQ,GetCoefficientWithAutomorphisms(x,d["RatsUsed"],d["AutUsed"]))
        ExplicitVector = QQ(ec["Coefficient"])*ExplicitVector
        ProducedIsotropicSequence.append(ExplicitVector.list())
    return ProducedIsotropicSequence

def GetDiagram(ec):
    # Given the initial data of curves and coefficients of an elliptic configuration, 
    # this function returns the Dynkin diagram type in the form [L, N], 
    # where L is a string with the letter identifying the Dynkin diagram,
    # and N is an integer identifying its length.
    # The program stops raising an exception in case the coefficients does not correspond to a valid Extended Dynkin diagram.
    # The program stops raising an exception if there are more than nine vertices, which is not allowed for an elliptic configuration on an Enriques surface.
    
    NumVertices=len(ec["EllipticConfiguration"])
    UniqueCoefficients=list(set([x[0] for x in ec["EllipticConfiguration"]]))
    UniqueCoefficients.sort()
    if (NumVertices > 9):
        raise Exception("Error: too many vertices in the Dynking diagram")
    if (UniqueCoefficients==[1]):
        return ["A",NumVertices-1]
    if (UniqueCoefficients==[1,2]):
        return ["D",NumVertices-1]
    if (UniqueCoefficients==[1,2,3]):
        return ["E",6]
    if (UniqueCoefficients==[1,2,3,4]):
        return ["E",7]
    if (UniqueCoefficients==[1,2,3,4,5,6]):
        return ["E",8]
    raise Exception("Error: the coefficients does not correspond to any valid Dynking diagram")
    
    
def extendedDynkinGraph(Type,Order):
    # extendedDynkinGraph generates the intersection matrix of an extended Dynkin diagram of Type "A", "D", or "E."
    # This function is the same as the one used in the CNDFinder code (github.com/rmoschetti/CNDFinder).
    # The result is returned as a dictionary with "Type", "Weights", and "Matrix" keys.
    # "Type" is a string identifying the type and order of the Dynkin diagram (e.g., "A4").
    # "Weights" is a list of weights for each node in the diagram.
    # "Matrix" is the intersection matrix for the diagram.
    # If Type is "A", it computes the matrix for type A Dynkin diagrams.
    # If Type is "D", it computes the matrix for type D Dynkin diagrams.
    # If Type is "E", it computes the matrix for type E Dynkin diagrams of orders 6, 7, or 8.

    Result={}
    Result["Type"]=Type + str(Order)
    
    if Type=="A":
        Result["Weights"]=[1]*(Order+1)
        
        if (Order==1):
            Result["Matrix"]=matrix([[-2,2],[2,-2]])
            
        else:
            Result["Matrix"]=(-2)*identity_matrix(Order+1)
            Result["Matrix"][0,-1]=1
            Result["Matrix"][-1,0]=1
            for i in range(0,Order):
                Result["Matrix"][i+1,i]=1
                Result["Matrix"][i,i+1]=1
        
    if Type=="D":
        Result["Weights"]=[2]*(Order+1)
        Result["Weights"][0]=1
        Result["Weights"][1]=1
        Result["Weights"][2]=1
        Result["Weights"][3]=1
        
        Result["Matrix"]=(-2)*identity_matrix(Order+1)
        Result["Matrix"][0,4]=1
        Result["Matrix"][1,4]=1
        Result["Matrix"][4,0]=1
        Result["Matrix"][4,1]=1
        Result["Matrix"][2,-1]=1
        Result["Matrix"][3,-1]=1
        Result["Matrix"][-1,2]=1
        Result["Matrix"][-1,3]=1
        for i in range(4,Order):
            Result["Matrix"][i+1,i]=1
            Result["Matrix"][i,i+1]=1
        
    if Type=="E":
        if (Order==6):
            Result["Weights"]=[1,2,3,2,1,2,1]
            Result["Matrix"]=matrix([
                            [-2, 1, 0, 0, 0, 0, 0],
                            [1, -2, 1, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 1, -2, 0, 0],
                            [0, 0, 1, 0, 0, -2, 1],
                            [0, 0, 0, 0, 0, 1, -2]
                        ])   
        if (Order==7):
            Result["Weights"]=[2,3,4,3,2,1,2,1]
            Result["Matrix"]=matrix([
                            [-2, 1, 0, 0, 0, 0, 0, 1],
                            [1, -2, 1, 0, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0, 0],
                            [0, 0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 0, 1, -2, 0, 0],
                            [0, 0, 1, 0, 0, 0, -2, 0],
                            [1, 0, 0, 0, 0, 0, 0, -2]
                        ])
        if (Order==8):
            Result["Weights"]=[2,4,6,5,4,3,2,3,1]
            Result["Matrix"]=matrix([
                            [-2, 1, 0, 0, 0, 0, 0, 0, 0],
                            [1, -2, 1, 0, 0, 0, 0, 0, 0],
                            [0, 1, -2, 1, 0, 0, 0, 1, 0],
                            [0, 0, 1, -2, 1, 0, 0, 0, 0],
                            [0, 0, 0, 1, -2, 1, 0, 0, 0],
                            [0, 0, 0, 0, 1, -2, 1, 0, 0],
                            [0, 0, 0, 0, 0, 1, -2, 0, 1],
                            [0, 0, 1, 0, 0, 0, 0, -2, 0],
                            [0, 0, 0, 0, 0, 0, 1, 0, -2]
                        ])
    return Result

def GetCoefficientWithAutomorphisms(CurveData,Rats,Aut):
    # In the data, a curve in an elliptic configuration is given by CurveData: an array of variable length such that
    # CurveData[0] is the coefficient of the curve in the elliptic configuration.
    # CurveData[1] is a placeholder for a length-10 vector in Rats.
    # CurveData[2, ...] is a placeholder for 10x10 matrices representing an automorphism in Aut.
    # For example, [2, "R14", "H0"] means 2*R14*H0.
    # This function returns explicitly the length-10 vector R14*H0 (coefficient excluded), 
    # which represents the curve R14*H0.
    
    if (len(CurveData)==2):
        return Rats[CurveData[1]]
    else:
        Ris=matrix(QQ,Rats[CurveData[1]])
        for A in CurveData[2::]:
            Ris = Ris * matrix(Aut[A])
        return Ris

def GetIntersection(A,B):
    # This function computes the intersection of two length-10 lists, A and B, interpreting them as vectors.
    # It uses the E10-matrix M as the intersection product.
    # The result is a scalar value representing the intersection of the two vectors.
    
    M=matrix([ [ -2, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], [ 0, -2, 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 1, -2, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 1, 0, 1, -2, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, -2, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, -2, 1, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 1, -2, 1, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 1, -2, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1, -2, 1 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 1, -2 ] ])
    
    return (matrix(QQ,A)*M*transpose(matrix(QQ,B)))[0,0]

def CheckEllipticConfigurationsInIsotropicSequence(d):
    # For all the elements in d["IsotropicSequence"] the function runs the following checks:
    # Check 1: Verify that the dual graph of the smooth rational curves is one of the allowed extended Dynkin diagrams: A1, ..., A8, D4, ..., D8, E6, E7, E8.
    # Check 2: Verify that the coefficients attached to each smooth rational curve match the coefficients of the associated extended Dynkin diagram.
    
    for Ell in d["IsotropicSequence"]:

        # Get the type of the diagram
        Type = GetDiagram(Ell)

        # Get the "correct" matrix for this type
        DynkDiagramType = extendedDynkinGraph(Type[0], Type[1])
        
        # Compute the intersection matrix of the curves in the elliptic configuration
        MatrixEllipticConfiguration = matrix(DynkDiagramType["Matrix"].ncols())
        for iv, v in enumerate(Ell["EllipticConfiguration"]):
            for iw, w in enumerate(Ell["EllipticConfiguration"]):
                MatrixEllipticConfiguration[iv, iw] = GetIntersection(
                        GetCoefficientWithAutomorphisms(v, d["RatsUsed"], d["AutUsed"]),
                        GetCoefficientWithAutomorphisms(w, d["RatsUsed"], d["AutUsed"])
                    )
        
        # Check 1: Compare the intersection matrices and verify that they are equal.
        # Notice that this is a sufficient condition; in principle, the matrices could be equal up to a permutation of rows and the corresponding columns.
        # In our case, there is no need to take into account such permutations because the smooth rational curves in the data we provided are sorted matching the ordering of the abstract extended Dynkin diagram.
        
        if (MatrixEllipticConfiguration != DynkDiagramType["Matrix"]):
            raise Exception("Error: one vertex has the wrong intersections to belong to the Dynkin diagram")
            
        CoefficientsEllipticConfiguration = [q[0] for q in Ell["EllipticConfiguration"]]
        
        # Check 2: Compare the coefficient lists and verify that they are equal. 
        # As for Check 1, this is also a sufficent condition.
        
        if (CoefficientsEllipticConfiguration != DynkDiagramType["Weights"]):
            raise Exception("The coefficient attached to the vertex of the Dynkin diagram is not correct")

def GenerateType(d):
    # This function generates a string representation of the isotropic sequence's type.
    # It combines the type identifiers for each elliptic configuration in the sequence
    # and includes "F" or "HF" based on the coefficients.
    
    Str = ""
    for _, x in enumerate(d["IsotropicSequence"]):
        Type = GetDiagram(x)
        Str = Str + Type[0] + str(Type[1])
        if (QQ(x["Coefficient"]) == 1/2):
            Str = Str + "F"
        else:
            Str = Str + "HF"
        Str = Str + ", "
    return Str[:-2]

    
def GenerateTex(d):
    # This function generates the LaTeX code to write the isotropic sequence.
    
    Str = ""
    Str = Str + "\\textbf{(" + str(d["no"]) + ")} "
    Str = Str + "$\\mathrm{nd}(Y_{" + str(d["no"]) + "})"
    
    if (d["cnd"] == 10):
        Str = Str + " = "
    else:
        Str = Str + " \\geq "
        
    Str = Str + str(d["cnd"]) + "$. "
    
    Str = Str + "Sequence: "
    
    for indiso, x in enumerate(d["IsotropicSequence"]):
        Str = Str + "$"
        
        if (QQ(x["Coefficient"]) == 1/2):
            Str = Str + "\\frac{1}{2}"

        Str = Str + "("
        for index, e in enumerate(x["EllipticConfiguration"]):
            if (e[0] != 1):
                Str = Str + str(e[0])
            Str = Str + "R_{"
            Str = Str + e[1][1::]
            Str = Str + "}"
            for indh, h in enumerate(e[2::]):
                Str = Str + " \\cdot "
                Str = Str + "H_{"
                Str = Str + h[1::]
                Str = Str + "}"
            if (index < len(x["EllipticConfiguration"]) - 1):
                Str = Str + "+"
        if (indiso < len(d["IsotropicSequence"]) - 1):
            Str = Str + ")$, "
        else:
            Str = Str + ")$.\\\\ \n"
        
    return Str

                
            
            
def CheckData(d):
    print("Checking", d["no"], end="")

    # Check that the coordinates of each isotropic vector in ["ProducedIsotropicSequence"] coincide
    # with the coordinates computed directly using the data in ["IsotropicSequence"]
    if (GetExplicitIsotropicSequence(d) != d["ProducedIsotropicSequence"]):
        raise Exception("The coefficients of the isotropic sequence already saved in the data DO NOT coincide with the coefficients computed directly using the original curves")

    # Check that the coordinates of each isotropic vector in ["ProducedIsotropicSequence"] are integers
    for v in d["ProducedIsotropicSequence"]:
        if (not all(isinstance(num, int) for num in v)):
            raise Exception("The coefficients of the isotropic sequence are not integers")

    # Check that the elliptic configurations with attached coefficient 1 cannot possibly be divided by 2
    for index, v in enumerate(d["IsotropicSequence"]):
        if (QQ(v["Coefficient"]) == 1 and all(num % 2 == 0 for num in d["ProducedIsotropicSequence"][index])):
            raise Exception("All the coordinates are even, and the elliptic configuration was not divided by two")

    # Check that the elliptic configurations claimed correspond to actual extended Dynkin diagrams
    CheckEllipticConfigurationsInIsotropicSequence(d)

    # Check that F_i \cdot F_j = 1 - \delta_{i,j}
    LengthSequence = len(d["ProducedIsotropicSequence"])
    DeltaMatrix = matrix([[1] * LengthSequence] * LengthSequence) - matrix.identity(LengthSequence)
    IntMatrix = matrix([[GetIntersection(A, B) for B in d["ProducedIsotropicSequence"]] for A in d["ProducedIsotropicSequence"]])
    if (DeltaMatrix != IntMatrix):
        raise Exception("Error: intersection of the elements of the isotropic sequences are not correct")

    # Check some additional data: the ["cnd"] should coincide with the length of the isotropic sequence
    if (d["cnd"] != LengthSequence):
        raise Exception("Error: cnd does not correspond")

    # Check some additional data: the ["type"] should coincide with the string returned by the function GenerateType
    if (GenerateType(d) != d["type"]):
        raise Exception("Error: the type is different than the one listed in the data")

    # Check some additional data: the ["tex"] should coincide with the string returned by the function GenerateTex
    if (GenerateTex(d) != d["tex"]):
        raise Exception("Error: the LaTeX is different than the one listed in the data")

    # If no exceptions had been raised, the data is valid, so print ok.
    print(" - ok")

    
import json
with open("Cnd_Data.txt") as json_file:
    Cnd_Data = json.load(json_file)


for d in Cnd_Data:
    CheckData(d)
 
     
