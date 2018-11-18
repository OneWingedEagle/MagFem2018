
package meshFactory;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.SourceDataLine;

import math.Complex;
import math.DFT;
import math.Vect;
import math.util;

public class pepe {

    public static void main(final String [] args) throws Exception {
    	
    	  File file=new File("C:/personal/Abdi/AbdiChacha17sec.wav");
    	  file=new File("C:/personal/Abdi/emergency003.wav");
    	  int totalFramesRead = 0;
    	  // somePathName is a pre-existing string whose value was
    	  // based on a user selection.
    	  Path path=file.toPath();
    	  
    	 
    	  
    	  byte[] data=Files.readAllBytes(path);
    	  int nbytes=1;
    	  int skip=1;
    	  Vect pp=new Vect(data.length/skip);
    	  int iz=0;
	      for(int i=0;i<pp.length;i+=skip){
	  
	    	  pp.el[iz++]=(double)(data[i] & 0xFF);
	      }
	      
	      double av=pp.sum()/pp.length;
	      
	      pp=pp.sub(new Vect().ones(pp.length).times(av));
	      util.plot(pp);
	      
/*	      Vect pp2=new Vect(pp.length);
	      int iz=0;
	      for(int i=0;i<pp2.length;i++){
	    	  pp2.el[i]=pp.el[i];
	      }
	   */
	      
	      Complex[] ff=DFT.fft(pp.el);
	      
/*	      Vect mag=new Vect(ff.length/10);
	      Vect f=new Vect(mag.length);
	      for(int i=1;i<mag.length;i++){
	    	  mag.el[i]=ff[i].norm();
	    	  f.el[i]=1.0*i/17;
	      }
	      util.plot(f,mag);*/
	      int kx=0;
	/*      for(int i=0;i<ff.length;i++){
	    	  if(1.0*i/17>100) break;
	    	  kx++;
	    	  ff[i]=ff[i].times(0) ;
	      }*/
  	    
	      for(int i=0;i<ff.length;i++){
	    	  double a=ff[i].norm()*(1+i*0e-3);
	    //	if(a>50000) a=50000;
	    	  double tt=ff[i].ang();
	    	  ff[i].re=a*Math.cos(tt);
	    	  ff[i].im=a*Math.sin(tt);
 	    	 
 	      }
	      
	      Vect mag=new Vect(ff.length);
	      Vect f=new Vect(mag.length);
	      for(int i=0;i<mag.length;i++){
	    	  mag.el[i]=ff[i+kx].norm();
	    	  f.el[i]=1.0*(kx+i)/17;
	      }
	      util.plot(f,mag);
	      
	      Complex[] ift=DFT.ifft(ff);
	      
	     Vect filtered=new Vect(pp.length);
	      for(int i=0;i<pp.length;i++){
	    	  filtered.el[i]=ift[i].re;
	      }
	    //  filtered.show();
		//  util.plot(filtered);
		 // v.show();
	      
	     // util.plot(pp2);
	   if(1.0>30)   
    	  try {
    	    AudioInputStream audioInputStream = 
    	      AudioSystem.getAudioInputStream(file);
    	    
    	    int numBytes = audioInputStream.available();
            System.out.println("numbytes: "+numBytes);
    	    int bytesPerFrame = 
    	      audioInputStream.getFormat().getFrameSize();
    	      if (bytesPerFrame == AudioSystem.NOT_SPECIFIED) {
    	      // some audio formats may have unspecified frame size
    	      // in that case we may read any amount of bytes
    	      bytesPerFrame = 2;
    	    } 
    	    // Set an arbitrary buffer size of 1024 frames.
    	    numBytes = 1024*16 * bytesPerFrame; 
    	    
    	 //   numBytes=1024*16*2;
    	    byte[] audioBytes = new byte[numBytes];
    	    
 
    	    Vect v=new Vect(audioBytes.length/bytesPerFrame);
    	    try {
    	      int numBytesRead = 0;
    	      int numFramesRead = 0;
    	      // Try to read numBytes bytes from the file.
    	      int iy=0;
    	      while ((numBytesRead = 
    	        audioInputStream.read(audioBytes)) != -1) {
    	        // Calculate the number of frames actually read.
    	        numFramesRead = numBytesRead / bytesPerFrame;
    	        totalFramesRead += numFramesRead;
    	        // Here, do something useful with the audio data that's 
    	        // now in the audioBytes array...
    	      }
  
    	    
   	      int ix=0;
    	      for(int i=0;i<v.length;i++){
    	
    		   
    	    	  
    	    	  v.el[i]=(double)(audioBytes[i] & 0xF);
    	      }
    	      
       	      util.plot(v);     
    		
       	   /*   
       	      Complex[] ff=DFT.fft(v.el);
    	    
    	      for(int i=0;i<ff.length;i++){
    	    	  double a=ff[i].norm()*(1+i*1e-3);
    	    	  double tt=ff[i].ang();
    	    	  ff[i].re=a*Math.cos(tt);
    	    	  ff[i].im=a*Math.sin(tt);
     	    	 
     	      }
    	      
    	      Complex[] ift=DFT.ifft(ff);
    	      
    	     Vect filtered=new Vect(v.length);
    	      for(int i=0;i<ff.length;i++){
    	    	  filtered.el[i]=ift[i].re;
    	      }
    	    //  filtered.show();
    		  util.plot(filtered);*/
    		 // v.show();
    	    } catch (Exception ex) { 
    	      // Handle the error...
    	    }
    	    
    	    
    
       	      
    	  } catch (Exception e) {
    	    // Handle the error...
    	  }
    	  
    	  
}    


} // End of the class //