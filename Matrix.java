package Linear_Algebra;

public class Matrix {
	
	private double[][] matrix;
	private double[][] point;
	
	public double[][] addition(double[][] m1, double[][] m2) {
		matrix = new double[m1.length][m2[0].length];
		if(m1.length*m1[0].length == m2.length*m2[0].length) {
			for(int i = 0; i < m1.length; i++) {
				for(int j = 0; j < m1[0].length;j++) {
					matrix[i][j] = m1[i][j] + m2[i][j];
				}
			}
			return matrix;
		} else {
			throw new IllegalArgumentException("Undefined Matrix");
		}
	}
	
	public double[][] subtraction(double[][] m1, double[][] m2) {
		matrix = new double[m1.length][m2[0].length];
		if(m1.length*m1[0].length == m2.length*m2[0].length) {
			for(int i = 0; i < m1.length; i++) {
				for(int j = 0; j < m1[0].length;j++) {
					matrix[i][j] = m1[i][j] - m2[i][j];
				}
			}
			return matrix;
		} else {
			throw new IllegalArgumentException("Undefined Matrix");
		}
	} 
	
	public double[][] scalarMultiplication(double[][] matrix, double scalar){
		this.matrix = new double[matrix.length][matrix[0].length];
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix[0].length; i++) {
				this.matrix[i][j] = matrix[i][j] * scalar; 
			}
		}
		return this.matrix;
	}
	
	   public boolean canMult(int rowLengthM1, int colLengthM2) {
		   if(rowLengthM1 == colLengthM2) {
			   return true;
		   } else {
			   return false;
		   }
	   }
	   
	   public double[][] hadamardProduct(double[][] m1, double[][] m2){
		   if(m1.length == m2.length && m1[0].length == m2[0].length) {
			   for(int i = 0; i < m1.length; i++) {
				   for(int j = 0; j < m1[0].length; j++) {
					   matrix[i][j] = m1[i][j] * m2[i][j];
				   }
			   }
			   return matrix;
		   } else {
			   throw new IllegalArgumentException("Matrices not same size");
		   }
	   }
	
	   public double[][] matrixMult(double[][] m1, double[][] m2) {
		   if(canMult(m1[0].length, m2.length) == false) {
			   throw new IllegalArgumentException("Does not Multiply");
		   } else {
			   matrix = new double[m1.length][m2[0].length];
			   double number = 0;
				   for(int i = 0; i < m1.length; i++) {
					   for(int j = 0; j < m2[0].length; j++) {
						   for(int k = 0; k < m1[0].length; k++) {
							   number += m1[i][k] * m2[k][j];
						   }
						   matrix[i][j] = number;
						   number = 0;
					   }
				   }
			   return matrix;
		   }
	   } 
	   
		public double[][] rotation2D(double angle) {
			matrix = new double[2][2];
			matrix[0][0] = Math.cos(angle);
			matrix[0][1] = -Math.sin(angle);
			matrix[1][0] = Math.sin(angle);
			matrix[1][1] = Math.cos(angle);
			return matrix;
		}
		
		/*
		 * matrices mult in order of Rx * Ry * Rz
		 */
		
		public double[][] rotation3DXYZ(double roll, double pitch, double yaw) {
			matrix = new double[3][3];
			matrix[0][0] = Math.cos(pitch)*Math.cos(yaw);
			matrix[0][1] = -Math.cos(pitch) * Math.sin(yaw);
			matrix[0][2] = Math.sin(pitch);
			matrix[1][0] = Math.sin(roll)*Math.sin(pitch)*Math.cos(yaw)+Math.cos(roll)*Math.sin(yaw);
			matrix[1][1] = -Math.sin(roll)*Math.sin(pitch)*Math.sin(yaw)+Math.cos(roll)*Math.cos(yaw);
			matrix[1][2] = -Math.sin(roll)*Math.cos(pitch);
			matrix[2][0] = -Math.cos(roll)*Math.sin(pitch)*Math.cos(yaw)+Math.sin(roll)*Math.sin(yaw);
			matrix[2][1] = Math.cos(roll)*Math.sin(pitch)*Math.sin(yaw)+Math.sin(roll)*Math.cos(yaw);
			matrix[2][2] = Math.cos(roll) * Math.cos(pitch);
			return matrix;
		}
		
		public void Rx(double theta) {
			matrix = new double[3][3];
			matrix[0][0] = 1;
			matrix[0][1] = 0;
			matrix[0][2] = 0;
			matrix[1][0] = 0;
			matrix[1][1] = Math.cos(theta);
			matrix[1][2] = -Math.sin(theta);
			matrix[2][0] = 0;
			matrix[2][1] = Math.sin(theta);
			matrix[2][2] = Math.cos(theta);
		}
		
		public void Ry(double phi) {
			matrix = new double[3][3];
			matrix[0][0] = Math.cos(phi);
			matrix[0][1] = 0;
			matrix[0][2] = Math.sin(phi);
			matrix[1][0] = 0;
			matrix[1][1] = 1;
			matrix[1][2] = 0;
			matrix[2][0] = -Math.sin(phi);
			matrix[2][1] = 0;
			matrix[2][2] = Math.cos(phi);
		}
		
		public void Rz(double psi) {
			matrix = new double[3][3];
			matrix[0][0] = Math.cos(psi);
			matrix[0][1] = -Math.sin(psi);
			matrix[0][2] = 0;
			matrix[1][0] = Math.sin(psi);
			matrix[1][1] = Math.cos(psi);
			matrix[1][2] = 0;
			matrix[2][0] = 0;
			matrix[2][1] = 0;
			matrix[2][2] = 1;
		}
		
		public double[][] stretch2D(double[][] point, double[][] stretchValues) {
			double[][] stretch = {{stretchValues[0][0],0},{0,stretchValues[1][0]}};
			this.point = new double[point.length][point[0].length];
			this.point = matrixMult(stretch, point); 
			return this.point;
		} //point and stretchValues must be horizontally represented
		
		public double[][] enlargement2D(double[][] point, double[][] enlargement) {
			double[][] Enlargement = {{enlargement[0][0],0},{0,enlargement[0][0]}};
			this.point = new double[point.length][point[0].length];
			this.point = matrixMult(Enlargement, point); 
			return this.point;
		}
		
		public double[][] stretch3D(double[][] point, double[][] stretchValues) {
			double[][] stretch = {{stretchValues[0][0],0,0},{0,stretchValues[1][0],0},{0,0,stretchValues[2][0]}};
			this.point = new double[point.length][point[0].length];
			this.point = matrixMult(stretch, point); 
			return this.point;
		} //point and stretchValues must be horizontally represented
		
		public double[][] enlargement3D(double[][] point, double[][] enlargement) {
			double[][] Enlargement = {{enlargement[0][0],0,0},{0,enlargement[0][0],0},{0,0,enlargement[0][0]}};
			this.point = new double[point.length][point[0].length];
			this.point = matrixMult(Enlargement, point); 
			return this.point;
		}
		
		public void printStructure(double[][] Matrix) {
			for(int i = 0; i < Matrix.length; i++) {
				for(int j = 0; j < Matrix[0].length; j++) {
					System.out.print(Matrix[i][j] + " ");
				}
				System.out.println();
			}
		}
		
		public double[][] identityMatrix(int size){
			matrix = new double[size][size];
			for(int i = 0; i < matrix.length; i++) {
				for(int j = 0; j < matrix[0].length; j++) {
					if(i == j) {
						matrix[i][j] = 1;
					} else {
						matrix[i][j] = 0;
					}
				}
			}
			return matrix;
		}
		
		public boolean isInverse(double[][] m1, double[][] m2) {
			matrix = new double[m1.length][m2[0].length];
			matrix = matrixMult(m1,m2);
			if(matrix.length != matrix[0].length) {
				throw new IllegalArgumentException("Matrix is not square, No inverse");
			} else {
				if(matrix == identityMatrix(matrix.length)) {
					return true;
				} else {
					return false;
				}
			}
		}
		
		public double[][] transpose(double[][] matrix){
			this.matrix = new double[matrix[0].length][matrix.length];
			for(int i = 0; i < matrix.length; i++) {
				for(int j = 0; j < matrix[0].length; j++) {
					this.matrix[j][i] = matrix[i][j];
				}
			}
			return this.matrix;
		}

}
