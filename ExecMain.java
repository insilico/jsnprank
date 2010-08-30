/*
 * Parse command-line options and run the SNPrank algorithm.
 * Authors:  Brett McKinney and Nick Davis
 * Email:  brett.mckinney@gmail.com, nick@nickdavis.name
 */
public class ExecMain {
	public static void main(String[] args) {
		
		if(args.length != 3) {
			System.out.println("Usage: java snprank matrix.txt gamma output.txt");
			System.exit(1);
		}
		
		ReadMatrix full_data = new ReadMatrix(args[0]);		
		full_data.snprank(full_data.getHeader(), full_data.getData(), args[1], args[2]);
	}
}
