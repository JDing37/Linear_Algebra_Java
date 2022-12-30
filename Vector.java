package Linear_Algebra;

public class Vector extends Matrix{
	
	public double magnitude;
	public double dotProduct;
	public double scalarProjection;
	public double[][] vector;
	public double inverseMagnitude;
	public double vectorProjection;
	public double vectorAngle;

	public double magnitude(double[][] vector) {
		magnitude  = 0;
		if(vector[0].length == 1) {
			for(int i = 0; i < vector.length; i++) {
				magnitude += Math.pow(vector[i][0],2);
			}
			magnitude = Math.sqrt(magnitude);
			return magnitude;
		} else {
			throw new IllegalArgumentException("Please enter vector as a column");
		}
	}
	
	public double dotProduct(double[][] v1, double[][] v2) {
		dotProduct = 0;
		if(v1[0].length == 1 && v2[0].length == 1 && v1.length == v2.length) {
			for(int i = 0; i < v1.length; i++) {
				dotProduct += (v1[i][0] * v2[i][0]);
			}
			return dotProduct;
		} else {
			throw new IllegalArgumentException("Please enter vectors as a column or enter vectors with same length");
		} 
	}
	
	public double scalarProjection(double[][] v1, double[][] v2) {
		magnitude = magnitude(v2);
		dotProduct = dotProduct(v1,v2);
		scalarProjection = dotProduct/magnitude;
		return scalarProjection;
	}
	
	public double[][] vectorProjection(double[][] v1, double[][] v2) {
		vector = new double[v1.length][1];
		magnitude = magnitude(v2);
		scalarProjection = scalarProjection(v1,v2);
		inverseMagnitude = 1 / magnitude;
		vectorProjection = scalarProjection * inverseMagnitude;
		for(int i = 0; i < v1.length; i++) {
			vector[i][0] = vectorProjection * v2[i][0];
		}
		return vector;
	}
	
	public double[][] unitVector(double[][] vector){
		vector = new double[vector.length][1];
		magnitude = magnitude(vector);
		inverseMagnitude = 1/magnitude;
		for(int i = 0; i < vector.length; i++) {
			this.vector[i][0] = inverseMagnitude * vector[i][0];
		}
		return this.vector;
	}
	
	public double angleBetweenVectors(double[][] v1, double[][] v2) {
		vectorAngle = Math.acos((dotProduct(v1,v2))/(magnitude(v1)*magnitude(v2)));
		return vectorAngle;
	}
	
	public double[][] crossProduct(double[][] v1, double[][] v2){
		vector = new double[v1.length][1];
		if(v1[0].length == 3 && v2[0].length == 3) {
			double[][] cutMatrix1 = {{v1[0][1],v1[0][2]},{v2[0][1],v2[0][2]}};
			double[][] cutMatrix2 = {{v1[0][0],v1[0][2]},{v2[0][0],v2[0][2]}};
			double[][] cutMatrix3 = {{v1[0][0],v1[0][1]},{v2[0][0],v2[0][1]}};
			vector[0][0] = determinant(cutMatrix1);
			vector[1][0] = -1*determinant(cutMatrix2);
			vector[2][0] = determinant(cutMatrix3);
			return vector;
		} else {
			throw new IllegalArgumentException("Please enter vectors as rows OR wrong dimension of vectors");
		}
	}
	
}
