from math import pi
def qpix(x,lambda_ = 1.547, Ldet = 2590. ):
    ''' DOCUMENT q(x)
   retuns q (in 1/A), x in mm
   
   q=2kSin(tth/2)=k*tth=2pi/lam*tth=(2 pi/lam) * (x/L)
   lambda=1.55 A
   L=2590mm have to look in the book for the exact one	
    '''

    lambda_=lambda_
    Ldet=Ldet
	
    qq=(2*pi/lambda_)*(x/Ldet);
    return qq	
 
