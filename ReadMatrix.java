import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

//this class to read a matrix file with variale names are in the first line.
//example matrix file:
// name1, name2, name3 ...
//
//


public class ReadMatrix {
	private String [] header;
	private double [][] data;
	
	public ReadMatrix(String file) {
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
//				data[index++] = Double.parseDouble(strArray);
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
		
		
	
	}//end of method readFile.
	
	
	public void pagerank(String[] name, double[][] data, String p, String outFile) {
		double [] c = sumCol(data);
		
		double [] r = sumRow(data);
		
		//find the indices of c array that element is not zero
		int [] k = findNonZeroIndex(c);
		
		double[][] D = sparse(k,k,reciprocal(c,k),data.length,data.length);
		double [][] e = getOnesMatrix(data.length, 1);
		
		double [][] I = speye(data.length,data.length);
		
		
		double [][] z_T = getMatrix(1, data.length, 1.0/data.length);
		for(int i=0; i<k.length; i++) {
			z_T[0][k[i]] = (1-Double.parseDouble(p))/data.length;
		}
		
		
		
		double[][] A = this.addMatrix(this.matrixMulti(this.matrixTimesDouble(data, Double.parseDouble(p)),D),
				this.matrixMulti(e, z_T));
		
		double [][] x = new double[e.length][e[0].length];
		for (int j = 0; j < e.length; j++) {
			x[j][0] = e[j][0]/data.length; 
		
		}
		
		double lambda=0.0;
		for(int i=0; i<5;i++) {
			lambda = 0.0;
			x=this.matrixMulti(A, x);
			
			
			for (int j = 0; j < x.length; j++) {
				lambda += x[j][0];
			}
			
			for (int j = 0; j < x.length; j++) {
				x[j][0] = x[j][0]/lambda;
			}
			
		}
		
		// otuput to file
        try {
			FileWriter fw = new FileWriter(outFile);
			BufferedWriter writer = new BufferedWriter(fw);
			writer.write("index, snp-rank, interaction in, interaction out, snp-name\n");
			for(int i=0; i<x.length; i++) {
				writer.write(i+1+", "+x[i][0]+", " + r[i] +", " +c[i]+", "+name[i]+"\n");
			}
			
			writer.close();
			fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		// Create a simple Bar chart
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (int j = 0; j < x.length; j++) {
        	dataset.setValue(x[j][0], "value", j+1+"");
        	System.out.println(x[j][0]);
		}
      
        ChartPanel myChartPanel = new ChartPanel(ChartFactory.createBarChart("SNP Rank",
              "SNP", "value", dataset, PlotOrientation.VERTICAL, false,
              true, false));
        
//      -- build a separate frame for the network show --
        JFrame f = new JFrame();
        f.setTitle("SNP Rank chart");
        f.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        f.setSize(800, 550);
        f.setLocationRelativeTo(null);
        f.getContentPane().setLayout(new BorderLayout());
        JPanel panel = new JPanel();
        panel.add(myChartPanel,BorderLayout.SOUTH);
       
        JScrollPane sp = new JScrollPane(panel);
        panel.setPreferredSize(new Dimension(780,550));
        panel.revalidate();
        f.getContentPane().add(sp);
        f.setVisible(true);
		//next 
		
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
		double [][] result = new double[a.length][a[0].length];
		for (int i = 0 ; i < result.length; i++ ){
			for(int j = 0; j < result[0].length; j++){
				   result[i][j] = a[i][j]*d;
			}
		}
		return result;
	}
	
	private double[][] matrixMulti(double[][] a, double [][] b){
		double [][] result = new double[a.length][b[0].length];
		for (int i = 0; i < a.length; i++) { //i is left Matrix's row
            for (int j = 0; j < b[0].length; j++) {//j is right Matrix's column
                for (int k = 0; k < a[0].length; k++) {//k is left Matrix's column
                    result[i][j] +=  a[i][k]* b[k][j];
                }//end k for loop :m1.column
            }//end j for loop :m2.column
        }//end i for loop :m1.row
		return result;
	}
	
	private double[][] speye(int m, int n){
		double[][] result = new double[m][n];
		//result[v1[i],v2[i]] = s[i], others is 0
		int index = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if(i==j) {
					result[i][j] = 1.0;
				}else {
					result[i][j] = 0.0;
				}
			}
		}
		return result;
	}
	
	/**
	 * 
	 * @param m : matrix rows
	 * @param n : matrix columns
	 * @return a m*n matrix filed with ones
	 */
	private double [][] getOnesMatrix(int m, int n) {
		double [][] result = new double[m][n];
		for (int i = 0; i < result.length; i++) {
			for (int j = 0; j < result[0].length; j++) {
				result[i][j] = 1.0;
			}
		}
		return result;
	}
	
	/**
	 * 
	 * @param m : matrix rows
	 * @param n : matrix columns
	 * @return a m*n matrix filed with v
	 */
	private double [][] getMatrix(int m, int n, double v) {
		double [][] result = new double[m][n];
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
		double [] result = new double[k.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = 1.0/c[k[i]];
		}
		
		return result;
	}
	
	private double [][] sparse(int[] v1, int[] v2, double[] s, int n, int m){
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
	
	
	
	private int [] findNonZeroIndex(double[] array) {
		ArrayList result = new ArrayList(array.length);
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
	
	private double [] sumCol(double[][] data) {
		double [] result = new double[data.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = 0.0;
			for (int j = 0; j < data[i].length; j++) {
				result[i] += (data[i][j]);
			}
		}
		return result;
	}
	
	private double [] sumRow(double[][] data) {
		double [] result = new double[data[0].length];
		for (int i = 0; i < result.length; i++) {
			result[i] = 0.0;
			for (int j = 0; j < data.length; j++) {
				result[i] += (data[j][i]);
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
