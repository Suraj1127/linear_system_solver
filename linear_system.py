#!/usr/bin/python3
import numpy as np

class LinearSystem:
    def __init__(self, no_of_equations, no_of_variables, A):
        """
        :param no_of_equations: number of equations of the linear system
        :param no_of_variables: number of variables of the linear system
        :A: matrix A containing coefficients of the linear system
        """
        self.no_of_equations = no_of_equations
        self.no_of_variables = no_of_variables
        self.A = A

    def row_echelonize(self):
        """
        reduce the matrix A into echelon form
        """
        last_pivot_row = -1
        for i in range(self.no_of_variables):
            if (self.A[last_pivot_row+1:,i] == np.zeros(self.no_of_equations-last_pivot_row-1)).all():
                continue
            pivot = self.A[last_pivot_row+1][i]
            if pivot == 0:
                for j in range(last_pivot_row+2,self.no_of_equations):
                    if self.A[j][i] != 0:
                        temp1 = np.copy(self.A[last_pivot_row+1])
                        temp2 = np.copy(self.A[j])
                        self.A[j], self.A[last_pivot_row+1] = temp1, temp2
                        break

            pivot = self.A[last_pivot_row+1][i]
            self.A[last_pivot_row+1] = self.A[last_pivot_row+1]/pivot
            last_pivot_row += 1

            for j in range(last_pivot_row+1, self.no_of_equations):
                if self.A[j][i] != 0:
                    self.A[j] = self.A[j] - self.A[j][i]*self.A[last_pivot_row]

    def reduced_row_echelonize(self):
        """
        reduce the row echelon form of A into reduced row echelon form
        """
        rank = 1
        pivot_columns = [0]
        for i in range(1, self.no_of_equations):
            try:
                pivot_column = list(self.A[i]).index(1)
            except ValueError:
                break
            rank += 1
            pivot_columns.append(pivot_column)
            for j in range(0, i):
                self.A[j] = self.A[j]-self.A[j,pivot_column]*self.A[i]
        self.rank = rank
        self.pivot_columns = pivot_columns
        self.free_columns = sorted(list(set(range(self.no_of_variables)) - set(pivot_columns)))

    def get_free_matrix(self):
        """
        returns the free matrix consisting of free variable columns of A
        """
        self.free_matrix = np.zeros((self.rank,self.no_of_variables-self.rank))
        for i in range(self.free_matrix.shape[1]):
            self.free_matrix[:,i] = self.A[:self.rank,self.free_columns[i]]
        return self.free_matrix

    def get_solution(self):
        """
        returns the column space matrix of the solution
        """
        self.get_free_matrix()
        if self.free_matrix.shape[1] == 0:
            return np.zeros((self.no_of_variables,1))
        self.solution_as_column_space = np.zeros((self.free_matrix.shape[0]+self.free_matrix.shape[1],len(self.free_columns)))
        identity_matrix = np.eye(len(self.free_columns))
        j = 0
        for i in self.pivot_columns:
            self.solution_as_column_space[i] = -self.free_matrix[j]
            j += 1
        j = 0
        for i in self.free_columns:
            self.solution_as_column_space[i] = identity_matrix[j]
            j += 1
        return self.solution_as_column_space

def main():
    # input for number of variables and equations
    print("Enter the number of variables:")
    no_of_variables = int(input())
    print("Enter the number of equations:")
    no_of_equations = int(input())

    # input for matrix A and setting to linear system
    A = np.zeros((no_of_equations, no_of_variables))
    print("Enter the coefficients of equations line by line separated by a single space:")
    print("\"\"\"\nFor example for 3 equations and 4 variables, the coefficient matrix A would be of the form:")
    print("4 4 2 3\n6 4 5 3\n3 6 4 2\n\"\"\"")
    for i in range(no_of_equations):
        A[i] = [int(number) for number in input().split()]

    linear_system = LinearSystem(no_of_equations, no_of_variables, A)
    linear_system.row_echelonize()
    linear_system.reduced_row_echelonize()
    print("\nSo, the column-space matrix of the solution is:\n{}".format(linear_system.get_solution()))

if __name__ == "__main__":
    main()


