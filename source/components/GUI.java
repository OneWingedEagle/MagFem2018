package components;


import java.awt.Color;
import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import main.Main;
import math.util;

public class GUI extends JFrame implements ActionListener{


	
	public		JPanel		panel,pp1;
	public JTextArea progressArea=new JTextArea(), paramArea=new JTextArea();
	public  TextField tfMeshFile,tfDataFile;
	public  TextField tfIterMax,tfErrorMax;
	public TextField[] tfX=new TextField[3];
	public Label[] lbX=new Label[3];
	public  Button Browse1,Browse2,bMainGUI,Run,bTerminate;
	public String dataFile,meshFile,fluxFilePath;
	
	public static int screenWidth,screenHeight;

		
	 
	 public GUI(String path) {

		
				Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
				screenWidth = (int)(screenSize.getWidth());
				screenHeight = (int)(screenSize.getHeight());
		
			
		 //======================  file paths	
				
				
				String meshFile= "";
				String dataFile= "";


			int	tag=22;	
				
				if(tag==0){
					//meshFile= System.getProperty("user.dir") + "\\quiz.txt";
				//	meshFile= System.getProperty("user.dir") + "\\quizFine.txt";
					//meshFile= System.getProperty("user.dir") + "\\condDular.txt";
				//	meshFile= System.getProperty("user.dir") + "\\NadaTest.txt";
					meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\model2D.txt";
				//	dataFile= System.getProperty("user.dir") + "\\dataQuiz.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataCondDular.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataCondSurf.txt";
					dataFile= System.getProperty("user.dir") + "\\inputs\\mesh\\dataModel2D.txt";
					
					//meshFile= System.getProperty("user.dir") + "\\bunLabyringth.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataCondDular.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataCondSurf.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataLabyringth.txt";
					//meshFile= System.getProperty("user.dir") + "\\prob20.txt";
					////dataFile= System.getProperty("user.dir") + "\\dataProb20.txt";
					
					//meshFile= System.getProperty("user.dir") + "\\coilcoreY1.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataCoilCoreY.txt";
					
					
					//meshFile= System.getProperty("user.dir") + "\\prismEl3layer.txt";
					//dataFile= System.getProperty("user.dir") + "\\dataPrismJ.txt";
					
					
					
					//meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\cube_unif.txt";
					//dataFile= System.getProperty("user.dir") + "\\inputs\\data\\data_unif.txt";
					//meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\model3D2.txt";
					//dataFile= System.getProperty("user.dir") + "\\inputs\\data\\dataModel3D2.txt";
				}
				
				else if(tag==1){
					
					meshFile= System.getProperty("user.dir") + "\\SIBC.txt";
					//meshFile= System.getProperty("user.dir") + "\\plungers\\plunger2Dv.txt";
/*					meshFile= System.getProperty("user.dir") + "\\plungers\\EMroger.txt";
						dataFile= System.getProperty("user.dir") + "\\dataEMroger.txt";*/
					dataFile= System.getProperty("user.dir") + "\\dataSIBC.txt";
						/*meshFile= System.getProperty("user.dir") + "\\bunMagFluid.txt";
						dataFile= System.getProperty("user.dir") + "\\dataMagFluid.txt";*/
						//dataFile= System.getProperty("user.dir") + "\\dataPlungerANS.txt";
				}
				
				else if(tag==-1){
					meshFile= System.getProperty("user.dir") + "\\stack.txt";
						dataFile= System.getProperty("user.dir") + "\\dataStack.txt";
						//dataFile= System.getProperty("user.dir") + "\\dataPlungerANS.txt";
						meshFile= System.getProperty("user.dir") + "\\cell.txt";
						dataFile= System.getProperty("user.dir") + "\\datamuMag.txt";
				}
				
				else if(tag==-2){
					meshFile= System.getProperty("user.dir") + "\\solidp.txt";
						dataFile= System.getProperty("user.dir") + "\\dataSolid.txt";
						//dataFile= System.getProperty("user.dir") + "\\dataPlungerANS.txt";
				}
				else if(tag==-3){
					meshFile= System.getProperty("user.dir") + "\\statStack.txt";
						dataFile= System.getProperty("user.dir") + "\\dataStatStack.txt";
						//dataFile= System.getProperty("user.dir") + "\\dataPlungerANS.txt";
				}
				else if(tag==-4){
					meshFile= System.getProperty("user.dir") + "\\statSolid.txt";
						dataFile= System.getProperty("user.dir") + "\\dataStatSolid.txt";
						//dataFile= System.getProperty("user.dir") + "\\dataPlungerANS.txt";
				}
				
				else if(tag==2){
					meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\mot4th2D.txt";
					dataFile= System.getProperty("user.dir") + "\\inputs\\data\\dataMot4th2D.txt";
				}
				
				else if(tag==3){
					meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\mot4th2DFine.txt";
					dataFile= System.getProperty("user.dir") + "\\inputs\\data\\dataMot4th2DFineJForFineMesh.txt";
					//	meshFile= System.getProperty("user.dir") + "\\caps\\motor8thFiner.txt";
					
					//dataFile= System.getProperty("user.dir") + "\\dataMechMotor8th.txt";
				}
				else if(tag==4){
					
					meshFile= System.getProperty("user.dir") + "\\statFrameSegPrismFinIns.txt";
					dataFile= System.getProperty("user.dir") + "\\dataMechStatFrameSeg3DLamin.txt";
				}
				
			else if(tag==5){
				meshFile= System.getProperty("user.dir") + "\\inputs\\mesh\\mot4th3D.txt";
				
				dataFile= System.getProperty("user.dir") + "\\inputs\\data\\dataMot4th3D.txt";
					
			/*	meshFile= System.getProperty("user.dir") + "\\bb.txt";
				
				dataFile= System.getProperty("user.dir") + "\\databb.txt";*/
				}
			else if(tag==6){

			//	meshFile= System.getProperty("user.dir") + "\\hexBeam.txt";
				meshFile= System.getProperty("user.dir") + "\\beam3D.txt";
				dataFile= System.getProperty("user.dir") + "\\dataBeam3D.txt";
			}
				
	else if(tag==7){
				
			meshFile= System.getProperty("user.dir") + "\\geoModel.txt";
				dataFile= System.getProperty("user.dir") + "\\dataGeo.txt";
			}
	else if(tag==8){
		
		meshFile= System.getProperty("user.dir") + "\\bento.txt";
	//	meshFile= System.getProperty("user.dir") + "\\bento3ang.txt";
			dataFile= System.getProperty("user.dir") + "\\dataBento.txt";
		}
				
	else if(tag==9){
		
		meshFile= System.getProperty("user.dir") + "\\plungerFull.txt";
			dataFile= System.getProperty("user.dir") + "\\dataPlungerFull.txt";
		}
				
	else if(tag==10){
		meshFile= System.getProperty("user.dir") + "\\plungers\\EMSlice.txt";
			dataFile= System.getProperty("user.dir") + "\\dataEM3D.txt";
	}
				
	else if(tag==11){
		
		meshFile= System.getProperty("user.dir") + "\\reactorLessAir.txt";
			dataFile= System.getProperty("user.dir") + "\\dataReactor.txt";
	
	
	}	
				
	else if(tag==12){
		
		meshFile= System.getProperty("user.dir") + "\\trans1f2D.txt";
			dataFile= System.getProperty("user.dir") + "\\dataTrans1f2DJ.txt";
			meshFile= System.getProperty("user.dir") + "\\bb.txt";
			dataFile= System.getProperty("user.dir") + "\\databb.txt";
		}
	else if(tag==13){
		
		meshFile= System.getProperty("user.dir") + "\\plungers\\plunger3D.txt";
		//meshFile= System.getProperty("user.dir") + "\\plungers\\react.txt";
	//	meshFile= System.getProperty("user.dir") + "\\hexaEl.txt";
		
			dataFile= System.getProperty("user.dir") + "\\dataPlungerHex.txt";
			//dataFile= System.getProperty("user.dir") + "\\dataReact.txt";
		}
				
	else if(tag==14){
		
		//meshFile= System.getProperty("user.dir") + "\\plungers\\plunger3Dsimp.txt";
		meshFile= System.getProperty("user.dir") + "\\plungers\\plunger3D0.txt";
	
			dataFile= System.getProperty("user.dir") + "\\dataPlunger3Dsimp3.txt";
		}
				
	else if(tag==15){
		
	//meshFile= System.getProperty("user.dir") + "\\quadEl.txt";
		meshFile= System.getProperty("user.dir") + "\\axi.txt";
	
			dataFile= System.getProperty("user.dir") + "\\dataAxi.txt";
		}
	else if(tag==16){
		
		//meshFile= System.getProperty("user.dir") + "\\plungers\\plunger3Dsimp.txt";
		meshFile= System.getProperty("user.dir") + "\\axi.txt";
	
			dataFile= System.getProperty("user.dir") + "\\dataAxihx.txt";
		}
	else if(tag==17){
		
	
		meshFile= System.getProperty("user.dir") + "\\caps\\motorFull.txt";
	
			dataFile= System.getProperty("user.dir") + "\\dataMechMotor.txt";
		}
	else if(tag==18){
		
		
		meshFile= System.getProperty("user.dir") + "\\caps\\motorHalf.txt";
	
			dataFile= System.getProperty("user.dir") + "\\dataMechMotorHalf.txt";
		}
	else if(tag==19){
		
		
		meshFile= System.getProperty("user.dir") + "\\caps\\motor8th.txt";
		//meshFile= System.getProperty("user.dir") + "\\caps\\motor8thFiner.txt";
			dataFile= System.getProperty("user.dir") + "\\dataMechMotor8th.txt";
			//dataFile= System.getProperty("user.dir") + "\\dataMot8th3D.txt";
		}
				
	else if(tag==20){
		
		
		meshFile= System.getProperty("user.dir") + "\\gears\\magGearRotOut.txt";
	
			dataFile= System.getProperty("user.dir") + "\\gears\\dataMagGear.txt";
			
	/*		meshFile= System.getProperty("user.dir") + "\\bunRing.txt";
			
			dataFile= System.getProperty("user.dir") + "\\dataRing.txt";*/
		}
	else if (tag==21){

		meshFile= "C:\\Users\\Me\\Desktop\\Gosh\\byMagFem\\bun2D.txt";
	
			dataFile= "C:\\Users\\Me\\Desktop\\Gosh\\byMagFem\\data2D.txt";
	}
	else if (tag==22){

		String folder="D:\\Works\\2018\\MagFemElectrost\\3D\\3coils\\roughmesh\\series-para";
	
		meshFile= folder + "\\bun.txt";
	
			dataFile=folder +  "\\data.txt";
	}
		
			
			//====================================
				

			panel = new JPanel(new FlowLayout(0,10,10));
			panel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			getContentPane().add(panel);
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			setTitle(" FEM Analysis : "+path);
		
			int width = (int)(.7*screenWidth);
			int height = (int)(.7*screenHeight);
			
			setSize(width,height);
			setLocation(10, 10);
	
		 
		    
	 //================================================================ redirecting console to text area
			
			int ww1=(int)(.65*width);
			int hh1=(int)(.7*height);
		
			progressArea.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			progressArea.setEditable(false);;
			progressArea.setFont(new Font("Arial", 0, (int)(10.0*screenWidth/1200)));
		
			progressArea.setBorder(BorderFactory.createLineBorder(Color.blue,1));
			   JScrollPane scrollPane = new JScrollPane(progressArea);
			   scrollPane.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createTitledBorder("Progress"),
						BorderFactory.createEmptyBorder(10,5,5,5)));
		  scrollPane.setPreferredSize(new Dimension(ww1,hh1));
		  
		  paramArea.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			paramArea.setEditable(false);;
			paramArea.setFont(new Font("Arial", 0, (int)(10.0*screenWidth/1200)));
			
			paramArea.setBorder(BorderFactory.createLineBorder(Color.blue,1));
			   JScrollPane scrollPane2 = new JScrollPane(paramArea);
			   scrollPane2.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createTitledBorder("Parameters"),
						BorderFactory.createEmptyBorder(10,5,5,5)));
		  scrollPane2.setPreferredSize(new Dimension((int)(.3*width),hh1));

	//=================================================================== configuring the panel		

				Label lbMeshFile=new  Label("Load Mesh"   , Label.RIGHT);
				Label lbDataFile=new  Label("Load Data"   , Label.RIGHT);
				
				int d1=(int)(80.*GUI.screenWidth/1200);
				int h1=(int)(20.*GUI.screenWidth/1200);
				lbMeshFile.setPreferredSize(new Dimension(d1,h1));
				lbDataFile.setPreferredSize(new Dimension(d1,h1));
				
				
				int d2=(int)(300.*GUI.screenWidth/1200);

				tfMeshFile=new TextField(meshFile);
				tfMeshFile.setPreferredSize(new Dimension(d2,h1));
							
				tfDataFile=new TextField(dataFile);
				tfDataFile.setPreferredSize(new Dimension(d2,h1));
				
				Browse1=new Button("Load Mesh");
				Browse1.setPreferredSize(new Dimension(d2/3,h1));
				Browse2=new Button("Load Data");
				Browse2.setPreferredSize(new Dimension(d2/3,h1));
				bMainGUI=new Button("Open Main GUI");
				bTerminate=new Button("Terminate");
				bTerminate.setPreferredSize(new Dimension(d2/3,h1));
				bMainGUI.setPreferredSize(new Dimension(d2,h1));
				
				Browse1.addActionListener(this);
				Browse2.addActionListener(this);
				bMainGUI.addActionListener(this);
				
				tfIterMax=new TextField("3000");
				tfIterMax.setPreferredSize(new Dimension(d2/5,h1));
				tfErrorMax=new TextField("1e-6");
				tfErrorMax.setPreferredSize(new Dimension(d2/5,h1));
				Label lbIterMax=new Label("ICCG Iteration max.");
				Label lbErrorMax=new Label(" Error max.");

				for(int i=0;i<3;i++){
					lbX[i]=new Label("");
					lbX[i].setPreferredSize(new Dimension(d2/5,h1));
				tfX[i]=new TextField("");
				tfX[i].setPreferredSize(new Dimension(d2/5,h1));
				}
			

				
				JPanel leftPanel = new JPanel(new GridLayout(6,1,10,10));
				leftPanel.setBorder(BorderFactory.createEmptyBorder(50,10,10,20));
				JPanel rightPanel = new JPanel(new FlowLayout(0,5,5));
				rightPanel.setBorder(BorderFactory.createEmptyBorder(50,5,5,0));
				
				panel.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
	
				Label lbsaveTo=new Label("Save output to :", Label.RIGHT);
				lbsaveTo.setFont(new Font("Arial", 1, 16));
				
				 Run=new Button("Run");
				Run.setBackground(Color.GREEN);
				Run.setPreferredSize(new Dimension(d2/3,h1));
				 
				Label empty1=new Label();
				empty1.setPreferredSize(new Dimension(50,3));
				Label empty2=new Label();
				empty2.setPreferredSize(new Dimension(50,3));

				JPanel filesPanel1 = new JPanel(new FlowLayout(0,10,10));
				
				filesPanel1.add(lbMeshFile);
				filesPanel1.add(tfMeshFile);
				filesPanel1.add(Browse1);
				filesPanel1.add(empty1);
				filesPanel1.add(Run);
				filesPanel1.add(bTerminate);
				
				JPanel filesPanel2 = new JPanel(new FlowLayout(0,10,10));
				filesPanel2.add(lbDataFile);
				filesPanel2.add(tfDataFile);
				filesPanel2.add(Browse2);
				filesPanel2.add(empty2);
			//	filesPanel2.add(bMainGUI);
				
				JPanel iterPanel = new JPanel(new FlowLayout(0,10,10));
				iterPanel.add(lbIterMax);
				iterPanel.add(tfIterMax);
				iterPanel.add(lbErrorMax);
				iterPanel.add(tfErrorMax);
				for(int i=0;i<3;i++){
				iterPanel.add(lbX[i]);	
				iterPanel.add(tfX[i]);	
				}
				
				
				JPanel filesPanel = new JPanel(new GridLayout(2,1,10,10));
				filesPanel.add(filesPanel1);
				filesPanel.add(filesPanel2);
				JPanel textPanel = new JPanel(new FlowLayout());
				textPanel.add(scrollPane);
				textPanel.add(scrollPane2);
				
				panel.add(filesPanel);
				panel.add(iterPanel);
				panel.add(textPanel);
				
				
			

	}
	 
		public static void main(String[] args){

			new Main();
		}
	 
	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==Browse1)
				getFile(1);
			else if(e.getSource()==Browse2)
				getFile(2);
			
	}
	 
		
		public void getFile(int i){
			FileDialog fd = new FileDialog(new Frame(),"Select bun  file",FileDialog.LOAD);
			fd.setVisible(true);
			fd.toFront();
			String Folder=fd.getDirectory();
			String File = fd.getFile();
			if(Folder!=null && File!=null)
			{
				if(i==1){
					meshFile=Folder+"\\"+File;
					tfMeshFile.setText(meshFile);
					fluxFilePath=Folder+"\\flux.txt";
				}
				else if(i==2){
					dataFile=Folder+"\\"+File;
					tfDataFile.setText(dataFile);
				}
		}
			fd.dispose();
		}
		
		public void writeLog(String logFilePath){
		//	String logFilePath = System.getProperty("user.dir") + "\\log.txt";
			
			try{
				String alltext=progressArea.getText();
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(logFilePath)));	
				for (String line : alltext.split("\\n")) 	
					pwBun.println(line);
				
				pwBun.close();
			}
			catch(IOException e){}
		}


}




