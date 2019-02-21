package math;


import java.awt.Color;

public class Mandelbrot {

    // return number of iterations to check if c = a + ib is in Mandelbrot set
    public static int mand(Complex z0, double c,int iterMax) {
        Complex z = z0;
        for (int t = 0; t < iterMax; t++) {
            if (z.norm() > c) return t;
            z = z.times(z).add(z0);
        }
        return iterMax;
    }

    public static void main(String[] args)  {
    	
    	args=new String[3];
    	args[0]=".0";
    	args[1]=".0";
    	args[2]="1";
        double xc   = Double.parseDouble(args[0]);
        double yc   = Double.parseDouble(args[1]);
        double size = Double.parseDouble(args[2]);

        int n   = 900;   // create n-by-n image
        int max = 255;   // maximum number of iterations
        
        double c=2;

        int iterMax=100;
                
        double x0 = -size;

        double y0 = -size;
        Complex z0 = new Complex(x0, y0);
        
        double dd=2*size/n;
         
        Picture picture = new Picture(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double dx = i*dd;
                double dy = j*dd;
                Complex dz = new Complex(dx, dy);
                Complex z = z0.add(dz.times(1));
                int val =iterMax- mand(z, c,iterMax);
     
                val=(int)(val*1./iterMax*max);
               // util.pr(val);
              //  if(gray>200 || gray<100) gray=max;
               // Color color = new Color(gray, gray, gray);
               Color color= getColor(val);
               
     
                picture.set(i, n-1-j, color);
            }
        }
        picture.show();
    }
    
    
    private static Color getColor(int c)
    {
        if(c<20) {
            int red = 0;
            int green = 0;
            int blue = 100 + (c * 155) / 20;
            return new Color(red, green, blue);
        }

        else if(c<25) 
        {
            int red = 0;
            int green = 0;
            int blue = 255;
            return new Color(red, green, blue);
        }

        else if(c<48) 
        {
            int red = 0;
            int green = ((c - 25) * 255) / 23;
            int blue = 255;
            return new Color(red, green, blue);
        }

        else if(c<73) 
        {
            int red = 0;
            int green = 190 + ((80 - c) * 65) / 32;
            int blue = ((73 - c) * 255) / 25;
            return new Color(red, green, blue);
        }

        else if(c<80) 
        {
            int red = 0;
            int green = 190 + ((80 - c) * 65) / 32;
            int blue = 0;
            return  new Color(red, green, blue);
        }

        else if(c<87) 
        {
            int red = 0;
            int green = 190 + ((c - 80) * 65) / 32;
            int blue = 0;
            return new Color(red, green, blue);
        }

        else if(c<112) 
        {
            int red = ((c - 87) * 255) / 25;
            int green = 190 + ((c - 80) * 65) / 32;
            int blue = 0;
            return  new Color(red, green, blue);
        }

        else if(c<160) 
        {
            int red = 255;
            int green = ((160 - c) * (160 - c) * 255) / 48 / 48;
            int blue = 0;
            return  new Color(red, green, blue);
        }
        
        return new Color(255,255,255);

}
}


