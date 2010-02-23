
public class ExecMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length != 3) {
			System.out.println("Usage: java ExecMain matrix_file P_value output_filename");
		}
		
		ReadMatrix myreadMatrix = new ReadMatrix(args[0]);
//		System.out.println("header: " );
//		for(int i=0;i<myreadMatrix.getHeader().length;i++) {
//			System.out.print(myreadMatrix.getHeader()[i]+" ,");
//		}
//		
//		System.out.println("\n\ndata: " );
//		for(int i=0;i<myreadMatrix.getData().length;i++) {
//			for (int j = 0; j < myreadMatrix.getData().length; j++) {
//				System.out.print(myreadMatrix.getData()[i][j]+" ,");
//			}
//			System.out.println("");
//		}
		
		myreadMatrix.pagerank(myreadMatrix.getHeader(), myreadMatrix.getData(), args[1], args[2]);

	}

}
