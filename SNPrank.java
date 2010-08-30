import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/*
 * Matrix methods and implementation of the SNPrank algorithm.
 * Authors:  Brett McKinney and Nick Davis
 * Email:  brett.mckinney@gmail.com, nick@nickdavis.name
 */
public class SNPrank {
	private String [] header;
	private double [][] data;
	
	public SNPrank(String file) {
		header = null;
		data = null;
		readFile(file);
	}
	
	private void readFile(String filename) {
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader br = new BufferedReader(fr);

			// strRow is used to read line from file
			String delimiter = "";
			String strRow = br.readLine();
			if(strRow.indexOf(',')>=0) {
				delimiter = "\\,";
			}else {
				delimiter = "\\s";
			}
			//set the header from first line of the input file
			this.setHeader(strRow.split(delimiter));
			
			data = new double[this.getHeader().length][this.getHeader().length];

			// set the data part from other lines of the input file
			int index = 0;
			while ((strRow = br.readLine()) != null) {
				if(strRow.indexOf(',')>=0) {
					delimiter = "\\,";
				}else {
					delimiter = "\\s";
				}
				String [] strArray = strRow.trim().split(delimiter);
				for(int i=0; i<data[index].length; i++) {
					data[index][i] = Double.parseDouble(strArray[i].trim());
				}
				index++;
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}
	
	
	public void snprank(String[] name, double[][] G, String gamma, String outFile) {
		// vector of column sums of G
		double [] colsum = sumCol(G);
		
		//find the indices of c array that element is not zero
		int[] colsum_nzidx = findNonZeroIndex(colsum);
		
		// sparse matrix where the nonzero indices are filled with 1/colsum[i]
		double[][] D = sparse(colsum_nzidx, colsum_nzidx, reciprocal(colsum, colsum_nzidx), G.length, G.length);
		double[][] e = getMatrix(G.length, 1, 1.0);
		
		// non-zero elements of colsum/d_j have (1 - gamma) in the numerator of 
		// the second term (Eq. 5 from SNPrank paper)
		double[][] T_nz = getMatrix(1, G.length, 1.0);
		for(int i=0; i <colsum_nzidx.length; i++) {
			T_nz[0][colsum_nzidx[i]] = 1 - Double.parseDouble(gamma);
		}
		
		double[][] T = addMatrix(matrixMulti(matrixTimesDouble(G, 
				Double.parseDouble(gamma)), D),	matrixMulti(e, T_nz));
		
		double[][] r = new double[e.length][e[0].length];
		for (int j = 0; j < e.length; j++) {
			r[j][0] = e[j][0]/G.length; 
		}
		
		double threshold = 1.0E-4;
		double lambda = 0.0;
		boolean converged = false;
		double[][] r_old = r;

		// if the absolute value of the difference between old and current r 
		// vector is < threshold, we have converged
		while(!converged) {
			lambda = 0.0; 
			r_old = r;
			r = matrixMulti(T, r);
			
			// sum of r elements
			for (int j = 0; j < r.length; j++) {
				lambda += r[j][0];
			}
			
			// normalize eigenvector r so sum(r) == 1
			for (int j = 0; j < r.length; j++) {
				r[j][0] = r[j][0]/lambda;
			}
			
			// check convergence
			for (int j = 0; j < r.length; j++){
				if (Math.abs(r[j][0] - r_old[j][0]) < threshold){
					converged = true;
				}
				
				else {
					converged = false;
					break;
				}
			}
		}
		
		// otuput to file, truncating values to 6 decimal places
        try {
			FileWriter fw = new FileWriter(outFile);
			BufferedWriter writer = new BufferedWriter(fw);
			writer.write("SNP\tSNPrank\tIG\n");
			for(int i=0; i<r.length; i++) {
				writer.write(name[i] + "\t" + String.format("%.6f", r[i][0]) + "\t" + String.format("%.6f", colsum[i]) + "\n");
			}
			
			writer.close();
			fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}		
	}
		
	private double[][] addMatrix(double[][] m1, double[][] m2){
		
		double[][] resultM = new double[m1.length][m1[0].length];
		for (int i = 0 ; i < resultM.length; i++ ){
			for(int j = 0; j < resultM[0].length; j++){
				   resultM[i][j] = m1[i][j] + m2[i][j];
			}
		}
			
		return resultM;
	}
	
	private double[][] matrixTimesDouble(double[][] a, double d){
		double[][] result = new double[a.length][a[0].length];
		for (int i = 0 ; i < result.length; i++ ){
			for(int j = 0; j < result[0].length; j++){
				   result[i][j] = a[i][j]*d;
			}
		}
		
		return result;
	}
	
	private double[][] matrixMulti(double[][] a, double [][] b){
		double[][] result = new double[a.length][b[0].length];
		for (int i = 0; i < a.length; i++) { //i is left Matrix's row
            for (int j = 0; j < b[0].length; j++) {//j is right Matrix's column
                for (int k = 0; k < a[0].length; k++) {//k is left Matrix's column
                    result[i][j] +=  a[i][k]* b[k][j];
                }//end k for loop :m1.column
            }//end j for loop :m2.column
        }//end i for loop :m1.row
		return result;
	}
	
	/**
	 * 
	 * @param m : matrix rows
	 * @param n : matrix columns
	 * @return a m*n matrix filed with v
	 */
	private double[][] getMatrix(int m, int n, double v) {
		double[][] result = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				result[i][j] = v;
			}
		}
		return result;
	}
	
	/**
	 * 
	 * @param k : a int array
	 * @return a double array filled with reciprocal of the int array 
	 */
	private double [] reciprocal(double [] c, int[] k) {
		double[] result = new double[k.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = 1.0/c[k[i]];
		}
		
		return result;
	}
	
	private double[][] sparse(int[] v1, int[] v2, double[] s, int n, int m){
		double[][] result = new double[n][m];
		//result[v1[i],v2[i]] = s[i], others is 0
		int index = 0;
		for (int i = 0; i < result.length; i++) {
			for (int j = 0; j < result[0].length; j++) {
				if(i==v1[index] && j==v2[index]) {
					result[i][j] = s[index++];
				}else {
					result[i][j] = 0;
				}
			}
		}
		
		return result;
	}
	
	private int[] findNonZeroIndex(double[] array) {
		ArrayList<Integer> result = new ArrayList<Integer>(array.length);
		for (int i = 0; i < array.length; i++) {
			if(array[i] != 0.0) {
				result.add(i);
			}
		}
		int[] tmp = new int[result.size()];
		for (int i = 0; i < tmp.length; i++) {
			tmp[i] = Integer.parseInt(result.get(i).toString());
		}
		return tmp;
	}
	
	private double[] sumCol(double[][] data) {
		double [] result = new double[data.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = 0.0;
			for (int j = 0; j < data[i].length; j++) {
				result[i] += (data[i][j]);
			}
		}
		return result;
	}
	
	/**
	 * @return the header
	 */
	public String[] getHeader() {
		return header;
	}
	/**
	 * @param header the header to set
	 */
	public void setHeader(String[] header) {
		this.header = header;
	}
	/**
	 * @return the data
	 */
	public double[][] getData() {
		return data;
	}
	/**
	 * @param data the data to set
	 */
	public void setData(double[][] data) {
		this.data = data;
	}
}
