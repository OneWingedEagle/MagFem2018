package math;

public class Penalty {

	
	public static void main(String[] args){
		

		int N=2;
		
		Mat M=new Mat(N,N);
		M.el[0][0]=2;
		M.el[0][1]=3;
		M.el[1][0]=3;
		M.el[1][1]=2;
		
		Vect b=new Vect(N);
		
		b.el[0]=4;
		b.el[1]=1;
		
		Vect cnst=new Vect(1,-1);
		
		
		Vect x2=lagrange(M,b,cnst);
	
		
		x2.hshow();
	}
	
	private static Vect lagrange(Mat M ,Vect b, Vect cns){
	
		
		int N=M.nRow;
				
				Mat M2=new Mat(N+1,N+1);
				M2.el[0][0]=2;
				M2.el[0][1]=3;
				M2.el[0][2]=cns.el[0];
				M2.el[1][0]=3;
				M2.el[1][1]=2;
				M2.el[1][2]=cns.el[1];;
				M2.el[2][0]=cns.el[0];;
				M2.el[2][1]=cns.el[1];;
				
			Vect b2=new Vect(N+1);
				for(int i=0;i<N;i++)
					b2.el[i]=b.el[i];

		//M2.show();
				Vect x2=M2.inv().mul(b2);

		return x2;

			}
		

	private static Vect penalty(Mat M ,Vect b, Vect cns){
	
		
		int N=M.nRow;

				
			Vect b2=new Vect(N+1);
				for(int i=0;i<N;i++)
					b2.el[i]=b.el[i];

		//M2.show();
				Vect x2=M.inv().mul(b2);

		return x2;

			}

}

