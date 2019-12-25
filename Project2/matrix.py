import math
from math import sqrt
import numbers

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """
        Calculates the determinant of a 1x1 or 2x2 matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")
        if self.h > 2:
            raise(NotImplementedError, "Calculating determinant not implemented for matrices largerer than 2x2.")
        v_ij = 1
        v_iNotj = 1
        if self.h ==1:
            return self.g[0]
        else:
            
            for i in range(self.h):
                for j in  range(self.w):
                    if i == j:
                        v_ij = v_ij*self[i][j]
                    else:
                        v_iNotj = v_iNotj*self[i][j]
            determ = v_ij - v_iNotj
            return determ
        

    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")
        sumT = 0
        for i in range(self.h):
            for j in range(self.w):
                if i == j:
                    sumT = sumT + self.g[i][j]
        return sumT

    def inverse(self):
        """
        Calculates the inverse of a 1x1 or 2x2 Matrix.
        """
        inv = []
        
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
        if self.h > 2:
            raise(NotImplementedError, "inversion not implemented for matrices larger than 2x2.")

        if self.h == 1:
            inv.append([1/self.g[0][0]])
        elif self.h == 2:
            if self.g[0][0] * self.g[1][1] != self.g[0][1] * self.g[1][0]:
                a = self.g[0][0]
                b = self.g[0][1]
                c = self.g[1][0]
                d = self.g[1][1]
                f = 1 / (a * d - b * c)
                a_new = a*f
                b_new = b*f*(-1)
                c_new = c*f*(-1)
                d_new = d*f
                inv = [ [d_new, b_new],[c_new, a_new] ]
            else:
                raise ValueError('inverting not possible.')
            
        return Matrix(inv)   

    def T(self):
        
        """
        Returns a transposed copy of this Matrix.
        """
        #transMatrix = zeroes(self.w, self.h)
        transMatrix = []
        for j in range(self.w):
            row_ma = []
            for i in range(self.h):
                row_ma.append(self.g[i][j]) 
            transMatrix.append(row_ma)    
        return Matrix(transMatrix)

    def is_square(self):
        return self.h == self.w

    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same") 
           
        matrixSum = []    
        for i in range(self.h):
            rowSum = []
            for j in range(self.w):
                value = self.g[i][j] + other.g[i][j]
                rowSum.append(value)
            matrixSum.append(rowSum)    
        return Matrix(matrixSum)

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        matrixNe = []
        for i in range(self.h):
            rowNe = []
            for j in range(self.w):
                 rowNe.append(self.g[i][j]*(-1))
            matrixNe.append(rowNe)        
        return Matrix(matrixNe)

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        matrixSub = []   
        for i in range(self.h):
            rowSub = []
            for j in range(self.w):
                value = self.g[i][j] - other.g[i][j]
                rowSub.append(value)
            matrixSub.append(rowSub)    
        return Matrix(matrixSub)

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        mul = []
        if self.w == other.h:
            transOther = other.T()
            for i in range(self.h):
                row_mul = []
                for j in range(other.w):
                    row_self = self.g[i]
                    row_transOther = transOther.g[j]
                    result = 0
                    for k in range(self.w):
                        result = result + (row_self[k] * row_transOther[k])
                    row_mul.append(result)    
                mul.append(row_mul)       
        else:
            raise(ValueError, "Can not be calculated")
        return Matrix(mul)    

    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        rmul = []
        if isinstance(other, numbers.Number):
            pass
            
            for i in range(self.h):
                row_rmul = []
                for j in range(self.w):
                    rmul_result = self.g[i][j]*other
                    row_rmul.append(rmul_result)
                rmul.append(row_rmul)    
        return Matrix(rmul)        
                
                