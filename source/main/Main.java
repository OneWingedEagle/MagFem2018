package main;
import math.*;

import static java.lang.Math.*;

import io.Console;

import java.awt.Color;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import javax.swing.JOptionPane;

import components.GUI;
import fem.Model;
import fem.RunCLN;
import fem.RunMag;
import fem.RunMagIPM;
import fem.RunMech;



	public class Main implements ActionListener, DropTargetListener{

	public GUI gui;

	private Model model;
	private Thread thread;
	private int nMesh=0,iterMax=2000;	
	private double  errMax;
	private String path = System.getProperty("user.dir");
	public  boolean console=true,dated=false;

	public Main()
	{		
		this.model=new Model();
		this.gui=new GUI(this.path);
		this.gui.Run.addActionListener(this);
		this.gui.bTerminate.addActionListener(this);
		this.gui.setVisible(true);
		//this.gui.Run.doClick();
	}	
	
	public static void main(String[] args){
	
		new Main();
	}

		
	public void loadSimple(){
		
		if(console){
		Console.redirectOutput(this.gui.progressArea);
		}
		this.model.meshFilePath=this.gui.tfMeshFile.getText();
		this.model.dataFilePath=this.gui.tfDataFile.getText();

		this.errMax=Double.parseDouble(this.gui.tfErrorMax.getText());
		this.iterMax=Integer.parseInt(this.gui.tfIterMax.getText());
		this.thread=new Thread(){
			long m1 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			@Override
			public void run(){

				double t_start= System.currentTimeMillis();
				
				Main.this.model.loadMesh(model.meshFilePath);

				prepare();

				
				model.loadData(model.dataFilePath);


	
				 if(model.magAnalysis && !model.mechAnalysis) {
					runMag(); 

					}
				else if(!model.magAnalysis && model.mechAnalysis) {

					runMech(); 

					}
				else if(!model.magAnalysis && !model.mechAnalysis) {

					runCLN(); 

					}

			double t_end= System.currentTimeMillis();
			System.out.format("Total cpu time (s): %10.1f\n",(t_end-t_start)/1000.);	

				String logFilePath = model.resultFolder+ "\\log.txt";
				gui.writeLog(logFilePath);
				gui.Run.setBackground(Color.green);
				gui.Run.setEnabled(true);

			}
		};
		this.thread.start();		
	}

	public void runMag(){

		//if(model.numberOfRegions==10 && model.motor){
			// RunMagGear mt=new RunMagGear();
			 //mt.runMag(model, this);
			//	}
	
		 if(model.numberOfRegions==17000 && model.motor){


		 RunMagIPM mt=new RunMagIPM();
		 mt.runMag(model, this);
		}else{
			 RunMag mt=new RunMag();

			 mt.runMag(model, this);
		}
		
}
	
	public void runCLN(){

		RunCLN cln=new RunCLN();

			 cln.run(model, this);
		
		
}
	
	public void runMech(){
		 RunMech mt=new RunMech();
		 mt.runMech(model, this);
		
}

		public void prepare(){
		Main.this.model.iterMax=Main.this.iterMax;
		Main.this.model.errCGmax=Main.this.errMax;

		
		
		try{

			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("").getAbsolutePath()+"\\_last_MagFEM_elected_path")));		
			pw.println(this.model.meshFilePath);
			pw.println(this.model.dataFilePath);

			pw.close();
		}
		catch(IOException e){}
		
		DateFormat dateFormat = new SimpleDateFormat("MM.dd.HH.mm.ss");
		Date date = new Date();
		String suff=dateFormat.format(date);

		
		this.model.resultFolder = new File(model.meshFilePath).getParentFile().getAbsolutePath();
	
		if(dated) this.model.resultFolder=this.model.resultFolder+suff;

		
	File folder = new File(this.model.resultFolder);
/*		if(folder.exists()){
			util.deleteDir(folder);
		}
		else
			folder.mkdir();
*/
		this.model.eddyFolder = 	this.model.resultFolder + "\\results";
		folder = new File( this.model.eddyFolder);
		if(folder.exists())
			util.deleteDir(folder);
		else
		folder.mkdir();

		this.model.fluxFilePath=System.getProperty("user.dir")+"\\flux.txt";
		this.model.eddyFilePath=System.getProperty("user.dir")+"\\eddy.txt";
		//this.gui.paramArea.setText("");
		if(console){
		Console.redirectOutput(this.gui.paramArea);
		Console.redirectOutput(this.gui.progressArea);
		}
		this.gui.Run.setBackground(Color.gray);
		this.gui.Run.setEnabled(false);
	}



	
	public void loadMesh(){	

				System.gc();
				long m1 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();

				Main.this.model.loadMesh(Main.this.model.filePath);

				long m2 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
				System.out.println("Used memory for  model setup: "+(m2-m1)/(1024*1024)+"MB");

					long m3 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
					System.out.println("Used memory for drawing mesh: "+(m3-m2)/(1024*1024)+"MB");
					System.out.println();
					System.out.println(" Number of regions: "+Main.this.model.numberOfRegions);
					System.out.println(" Total number of elements: "+Main.this.model.numberOfElements);
					System.out.println(" Total number of nodes: "+Main.this.model.numberOfNodes);

					System.out.println();

	}		


	public void loadFlux(){

		boolean result=this.model.loadFlux( this.model.filePath);

		}	


	public void loadPotential(){
		
	

		boolean result=this.model.loadPotential( this.model.filePath);
	

		if(result) this.model.setB();


	}	

	public void loadStress(){


		boolean result=this.model.loadStress( this.model.filePath);
	
	}	

	public void loadNodalField(int mode){

		boolean result=this.model.loadNodalField( this.model.filePath, mode);

	}	

	public void paintDisplacement(){

		boolean result=this.model.loadNodalField( this.model.filePath, 0);
			
	}	


	public String readFirst(String filePath){
		try{

			Scanner scr=new Scanner(new FileReader(filePath));
			String first= scr.next();
			scr.close();
			return first;
		}
		catch(IOException e){return null;}

	}


	

	@Override
	public void actionPerformed(ActionEvent e)
	{	

		 if(e.getSource()==this.gui.Run){
			this.gui.Run.setBackground(Color.gray);
			this.gui.Run.setEnabled(false);

			loadSimple();
			//loadMotor();
		}
		else if(e.getSource()==this.gui.bTerminate){
			this.model.solver.terminate(true);

		}

	
	}


	@Override
	public void dragEnter(DropTargetDragEvent dtde) 
	{}

	@Override
	public void dragExit(DropTargetEvent dte) 
	{ }

	@Override
	public void dragOver(DropTargetDragEvent dtde) 
	{ }

	@Override
	public void dropActionChanged(DropTargetDragEvent dtde) 
	{  }

	@Override
	public void drop(DropTargetDropEvent dtde) {	   
		try {
			Transferable tr = dtde.getTransferable();
			DataFlavor[] flavors = tr.getTransferDataFlavors();
			
			for (int i = 0; i < flavors.length; i++) {
				if (flavors[i].isFlavorJavaFileListType()) {
					dtde.acceptDrop(DnDConstants.ACTION_COPY);
					List list = (List) tr.getTransferData(flavors[i]);
					this.model.filePath=list.get(i).toString();
				
					System.out.println("Dropped File: "+ this.model.filePath);
					String str=readFirst( this.model.filePath);

					if(str.equals("hexahedron") ||
							str.equals("triangle")||
							str.equals("quadrangle")||
							str.equals("prism")){
						this.thread=new Thread(){
							@Override
							public void run(){
								loadMesh(); 
							}
						};
						this.thread.start();
						this.nMesh++;
					}
					else if (str.equals("flux")){
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									loadFlux(); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}
					else if (str.equals("vPot")){
						
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									loadPotential(); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}

					else if (str.equals("stress")){
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									loadStress(); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}

					else if (str.equals("displacement")){
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									//paintDisplacement(); 
									loadNodalField(-1); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}


					else if (str.equals("force_reluc")){
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									loadNodalField(1); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}

					else if (str.equals("force_ms")){
						if(this.nMesh>0){
							this.thread=new Thread(){
								@Override
								public void run(){
									loadNodalField(2); 
								}
							};

							this.thread.start();

						}
						else{
							String msg="No mesh loded yet.";
							JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
						}
					}
					
					else {
						String msg="Invalid input file.";
						JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
					}

					
				}
				dtde.dropComplete(true);
			
			}

			System.out.println("Drop failed: " + dtde);
			dtde.rejectDrop();
		} catch (Exception e) {
			e.printStackTrace();
			dtde.rejectDrop();
		}
	}

	public void wait(int ms){
		try {
			Thread.sleep(ms);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public double getDiffMax(Vect[] u, Vect[] v){
		double diff;
		double diffMax=0;
		for(int i=0;i<u.length;i++){
			diff=u[i].sub(v[i]).norm();
			if(diff>diffMax)
				diffMax=diff;
		}

		return diffMax;
	}

	public double getErrorMax(Vect[] u, Vect[] v){
		double diff;
		double diffMax=0;
		double vmax=0;
		for(int i=0;i<u.length;i++){
			double vn=v[i].norm();
			if(vn>vmax) vmax=vn;
			diff=u[i].sub(v[i]).norm();
			if(diff>diffMax)
				diffMax=diff;
		}

		if(vmax==0)
			return 0;

		return diffMax/vmax;
	}





}

