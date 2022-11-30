package Linear_Algebra;

public class Vector extends Matrix{
	
	public double magnitude;
	public double dotProduct;
	public double scalarProjection;
	public double[][] vector;
	public double inverseMagnitude;
	public double vectorProjection;

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
		magnitude = magnitude(vector);
		inverseMagnitude = 1/magnitude;
		for(int i = 0; i < vector.length; i++) {
			this.vector[i][0] = inverseMagnitude * vector[i][0];
		}
		return this.vector;
	}
	
}
