package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Scanner;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class EMSolOutputReader {

	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		EMSolOutputReader x=new EMSolOutputReader();

		//Vect  NRError=x.loadOutput();
		
		
		String file="C:\\Works\\Problems and Models\\Large-scale Test\\FieldSources\\transient_linear_series _current_input\\MPI para 8\\refine 0\\output";
		
		//file="C:\\Users\\Hassan Ebrahimi\\Desktop\\new modules\\new\\output";


		Vect  i=x.loadSourceCurrent(file,1);
		Vect  v=x.loadSourceVoltage(file,1);
	
		//i.times(1).show();
		//util.plot(v);
		//util.plot(i);
	v.show();
		
	//	Vect y=x.loadMagEnergy(file,6); y.show();
		//Vect y=x.loadMagFlux(file,1);y.times(1000).show();
	

		//Vect  NRError=x.loadEnergy();
	}

	public Vect loadOutput(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_iso\\output";
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\output";
	//	String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\output";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		
		Vect tempICCGEr=new Vect(100000);
		Vect tempNREr=new Vect(100000);
		Vect tempdB=new Vect(100000);
		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** ICCG converges ***")){}
			if(line==null) break;
			line=br.readLine();
			sp=line.split(regex);
			tempICCGEr.el[ix]=Double.parseDouble(sp[sp.length-1]);
			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex);
			tempNREr.el[ix]=Double.parseDouble(sp[sp.length-1]);
			
			line=br.readLine();

			sp=line.split(regex);
			tempdB.el[ix]=Double.parseDouble(sp[4]);
			
			ix++;

			
			}
	br.close();
	fr.close();
	
	Vect errs=new Vect(ix);
	Vect NRerr=new Vect(ix);
	Vect dBs=new Vect(ix);
	for(int i=0;i<ix;i++){
		errs.el[i]=tempICCGEr.el[i];
		NRerr.el[i]=tempNREr.el[i];
		dBs.el[i]=tempdB.el[i];
	}
	
	util.plot(dBs);
	//util.plot(NRerr);
	//errs.show();
	return errs;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadCheck(){


		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_iso\\output";
		//String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Large model_Angs\\output";
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\check";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		
		Vect tempICCGEr=new Vect(100000);
	
		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** ICCG iteration")){}
			if(line==null) break;
			line=br.readLine();
		
			while(!(line=br.readLine()).startsWith("*")){
		
				sp=line.split(regex);
				//util.show(sp);
				tempICCGEr.el[ix]=Double.parseDouble(sp[1]);
				ix++;
			}
			
			
		
			
			continue;
			
			}
	br.close();
	fr.close();
	
	Vect errs=new Vect(ix);
	//Vect NRerr=new Vect(ix);
	//Vect dBs=new Vect(ix);
	for(int i=0;i<ix;i++){
		errs.el[i]=tempICCGEr.el[i];
		//NRerr.el[i]=tempNREr.el[i];
		//dBs.el[i]=tempdB.el[i];
	}
	
	util.plot(errs);
	//util.plot(NRerr);
	//errs.show();
	return errs;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	
	}
	
	public Vect loadEnergy(){

	
		String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\ringCompositAngDep\\output";

		
		Vect tempEnrgy=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic energy ")){}
			if(line==null) break;
			line=br.readLine();
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			tempEnrgy.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	Vect energy=new Vect(ix);

	for(int i=0;i<ix;i++){
		energy.el[i]=tempEnrgy.el[i];

	}
	
	util.plot(energy);
	//util.plot(NRerr);
	energy.times(1e6).show();
	return energy;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadSourceCurrent(String file, int sourceID){

	
		
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources")){}
			if(line==null) break;
			line=br.readLine();
			for(int i=0;i<sourceID;i++)
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	Vect current=new Vect(ix);

	for(int i=0;i<ix;i++){
		current.el[i]=powerSourceCurrent.el[i];

	}
	
	//util.plot(current);
	//util.plot(NRerr);

	return current;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	
	public Vect loadSourceVoltage(String file,int sourceID){

	
		Vect powerSourceCurrent=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources")){}
			if(line==null) break;
			line=br.readLine();
			for(int i=0;i<sourceID;i++)
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			powerSourceCurrent.el[ix]=Double.parseDouble(sp[3]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	Vect voltage=new Vect(ix);

	for(int i=0;i<ix;i++){
		voltage.el[i]=powerSourceCurrent.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return voltage;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	
	
	public Vect loadMagFlux(String file,int loopID){

		
		
		Vect magFluxTemp=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*            Magnetic fluxes of flux loops")){}
			if(line==null) break;
			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
		
			for(int i=0;i<loopID;i++)
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			magFluxTemp.el[ix]=Double.parseDouble(sp[2]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	Vect magFlux=new Vect(ix);

	for(int i=0;i<ix;i++){
		magFlux.el[i]=magFluxTemp.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return magFlux;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	

	
public Vect loadMagEnergy(String file,int loopID){

		
		
		Vect magEnergyTemp=new Vect(10000);

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic energy")){}
			if(line==null) break;
			line=br.readLine();
			for(int i=0;i<loopID;i++)
			line=br.readLine();
		//	line=br.readLine();
			sp=line.split(regex);
			magEnergyTemp.el[ix]=Double.parseDouble(sp[3]);
			
			ix++;

			
			}
			
	br.close();
	fr.close();
	
	Vect magFnergy=new Vect(ix);

	for(int i=0;i<ix;i++){
		magFnergy.el[i]=magEnergyTemp.el[i];

	}
	
	//util.plot(voltage);
	//util.plot(NRerr);

	return magFnergy;
		}

		catch(IOException e){System.err.println("Error in loading output file.");
		return null;
		}

		
	}	


}
