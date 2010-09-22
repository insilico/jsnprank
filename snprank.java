/*
 * SNPrank - single nucleotide polymorphism (SNP) ranking algorithm
 *
 * Uses a GAIN file, together with a damping factor gamma, 
 * (default is .85), to compute SNPrank scores.
 * Prints a series of rows containing the SNP name, SNPrank score,
 * information gain, sorted in descending order by SNPrank.
 * 
 * Authors:  Brett McKinney <brett.mckinney@gmail.com>
             Nick Davis <nick@nickdavis.name>
 */
public class snprank {
	public static void main(String[] args) {
		
		if(args.length != 3) {
			System.out.println("Usage: java snprank matrix.txt gamma output.txt");
			System.exit(1);
		}
		
		SNPrankData full_data = new SNPrankData(args[0]);		
		full_data.snprank(full_data.getHeader(), full_data.getData(), args[1], args[2]);
	}
}
