class LocalisationRender{
	int frame;
	int origX;
	int origY;
	float error;
	float noise;
	float bkg; //background
	float signal;
	float angle;
	float x;
	float y;
	float xsd;
	float ysd;
	float precision;
	float delta;
	float a; //fit parameter see Gaussian equation
	float b;
	float c;
	float d;
	float r2;
	boolean use; 
	float z0z1;
	float lambdaRef;
	float lambda;

	
	  public LocalisationRender(int frame, int origX, int origY, float error, float noise, float bkg, float signal, float angle, float x, float y, float xsd, float ysd, float precision, float delta){
		    this.frame = frame;
	    	this.origX = origX;
	    	this.origY = origY;
	    	this.error = error;
	    	this.noise = noise;
	    	this.bkg = bkg;
	    	this.signal = signal;
	    	this.angle = angle;
	    	this.x = x;
	    	this.y = y;
	    	this.xsd = xsd;
	    	this.ysd = ysd;
	    	this.precision = precision;
	    	this.delta = delta;
	  }
	
    public LocalisationRender(int frame, int origX, int origY,float error, float noise, float bkg, float signal, float angle, float x, float y, float xsd, float ysd, float precision, float delta, float lambdaRef){
    	//localizations from the calibration file
    	this.frame = frame;
    	this.origX = origX;
    	this.origY = origY;
    	this.error = error;
    	this.noise = noise;
    	this.bkg = bkg;
    	this.signal = signal;
    	this.angle = angle;
    	this.x = x;
    	this.y = y;
    	this.xsd = xsd;
    	this.ysd = ysd;
    	this.precision =precision;	
    	this.delta = delta;
    	this.lambdaRef = lambdaRef;
    	
    }
	
    public LocalisationRender(int frame, int origX, int origY,float error, float noise, float bkg, float signal, float angle, float x, float y, float xsd, float ysd, float precision, float delta, float lambda, boolean use){
    	this.frame = frame;
    	this.origX = origX;
    	this.origY = origY;
    	this.error = error;
    	this.noise = noise;
    	this.bkg = bkg;
    	this.signal = signal;
    	this.angle = angle;
    	this.x = x;
    	this.y = y;
    	this.xsd = xsd;
    	this.ysd = ysd;
    	this.precision =precision;	
    	this.delta = delta;
    	this.lambda = lambda;
    	this.use = use;
    }
    
    public void use(){
    	if( a!=0){
    	boolean condAmpl = b > 0 && a > 0;
    	boolean condCentre = 0 < c && c < 20;
    	boolean condWidth = 0 < d && d < 10;
    	this.use = condAmpl && condCentre && condWidth;
    	}else{
        this.use = false;
    	}
    }
    
    
    //getters 
    public int getFrame(){
    	return(frame);
    }
    public float getNoise(){
    	return(noise);
    }
    public float getBkg(){
    	return(bkg);
    }
    public float getSignal(){
    	return(signal);
    }
    public float getAngle(){
    	return(angle);
    }
    public float getX(){
    	return(x);
    }
    public float getY(){
    	return(y);
    }
    public float getXsd(){
    	return(xsd);
    }
    public float getYsd(){
    	return(ysd);
    }
    public float getDelta(){
    	return(delta);
    }
    public float getPrecision(){
    	return(precision);
    }
    public float getLambda(){
    	return(lambda);
    }
    public float getA(){
    	return(a);
    }
    public float getB(){
    	return(b);
    }
    public float getC(){
    	return(c);
    }
    public float getD(){
    	return(d);
    }
    public float getR2(){
    	return(r2);
    }
    public boolean getUse(){
    	return(use);
    }
    public float getZ0z1(){
    	return(z0z1);
    }
    public float getLambdaRef(){
    	return(lambdaRef);
    }
    //other setters???
    public void setA(float x){
    	this.a = x;
    }
    public void setB(float x){
    	this.b = x;
    }
    public void setC(float x){
    	this.c = x;
    }
    public void setD(float x){
    	this.d = x;
    }
    public void setR2(float x){
    	this.r2 = x;
    }
    public void setZ0z1(float x){
    	this.z0z1 = x;
    }
    public void setLambda(float x){
    	this.lambda= x;
    }
    public void setUseIt(boolean a){
    	this.use = a;
    }
    
    public String println(){
    	String a = ""+frame+"\t"+origX+"\t"+ origY +"\t"+ error +"\t"+ noise +"\t"+ bkg +"\t"+ signal +"\t"+angle+"\t"+ x +"\t"+ y +"\t"+ xsd +"\t"+ ysd +"\t"+ precision +"\t"+ z0z1 +"\t"+ lambda +"\t"+ use + "\n";
    	return(a);
    }
    
}//fin de la classe localisation

/*
 * @ Frame(0)	origX(1)	origY(2)	origValue(3)	Error(4)	Noise(5)
 * @ Background(6)	Signal(7)	Angle(8)	X(9)	Y(10)	X SD(11)	Y SD(12)	Precision(13)
 * Z0Z1(14) lamba_nm(15)
 */
