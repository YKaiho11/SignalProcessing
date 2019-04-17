#pragma once

class comp{
public:
    double re,im;
    
    comp(double r,double i){
        re=r;im=i;
    }
    comp(){
        re=0;im=0;
    }
    
    float real()const{return re;}
    float imag()const{return im;}
    
    void show(){
        printf("%f + %fi\n",re,im);
    }
    
    
    comp& operator=(comp& z){
        this->re=z.re;        this->im=z.im;
        return *this;
    }
    
    comp& operator+=(comp& z){
        this->re+=z.re;this->im+=z.im;
        return *this;
    }
    
    comp& operator+=(float x){
        this->re+=x;
        return *this;
    }
    
    comp& operator-=(comp& z){
        this->re-=z.re;this->im-=z.im;
        return *this;
    }
    
    comp& operator-=(float x){
        this->re-=x;
        return *this;
    }
    
    comp& operator*=(comp& z){
        float re_tmp = this->re*z.re - this->im*z.im;
        float im_tmp = this->re*z.im + this->im*z.re;
        this->re = re_tmp;    this->im = im_tmp;
        return *this;
    }
    
    comp& operator*=(float x){
        this->re *=x;    this->im *=x;
        return *this;
    }
    
    comp& operator/=(comp& z){
        float norm_z = z.re*z.re+z.im*z.im;
        float re_tmp = (this->re * z.re + this->im * z.im) / norm_z;
        float im_tmp = (this->im * z.re - this->re * z.im) / norm_z;
        this->re = re_tmp;    this->im = im_tmp;
        return *this;
    }
    
    comp& operator/=(float x){
        this->re /= x;    this->im /= x;
        return *this;
    }
};

//Complex+Complex
inline comp operator+(const comp& x,const comp& y){
    return comp(x.real() + y.real(),x.imag() + y.imag());
}
//Complex+real
inline comp operator+(const comp& x,const float y){
    return comp(x.real() + y, x.imag());
}
//real+Complex
inline comp operator+(const float x,const comp& y){
    return comp(x + y.real(), y.imag());
}
//Complex-Complex
inline comp operator-(const comp& x,const comp& y){
    return comp(x.real() - y.real(), x.imag() - y.imag());
}
//Complex-real
inline comp operator-(const comp& x,const float y){
    return comp(x.real() - y, x.imag());
}
//real-Complex
inline comp operator-(const float x,const comp& y){
    return comp(x - y.real(), -y.imag());
}
//Complex*Complex
inline comp operator*(const comp& x,const comp& y){
    return comp (
                    x.real() * y.real() - x.imag() * y.imag()    ,
                    x.real() * y.imag() + x.imag() * y.real()    );
}
//Complex*float
inline comp operator*(const comp& x,float y){
    return comp (
                    x.real() * y    ,
                    x.imag() * y    );
}
//float*Complex
inline comp operator*(float x,const comp& y){
    return comp (
                    x * y.real()    ,
                    x * y.imag()    );
}
//Complex/Complex
inline comp operator/(const comp& x,const comp& y){
    float norm_y = y.real()*y.real()+y.imag()*y.imag();
    return comp (
                    (x.real() * y.real() + x.imag() * y.imag()) / norm_y    ,
                    (x.imag() * y.real() - x.real() * y.imag()) / norm_y    );
}
//Complex/float
inline comp operator/(const comp& x,float y){
    return comp (
                    x.real() / y    ,
                    x.imag() / y    );
}
//float/Complex
inline comp operator/(float x,const comp& y){
    float norm_y = y.real() * y.real() + y.imag() * y.imag();
    return comp (
                    x * y.real() / norm_y        ,
                    -x * y.imag() / norm_y    );
}

inline bool operator==(const comp& x,const comp& y){
    return ( x.real() == y.real() && x.imag() == y.imag() );
}
inline bool operator==(const comp& x,const float& y){
    return ( x.real() == y && x.imag() == 0 );
}
inline bool operator==(const float& x,const comp& y){
    return ( x == y.real() && y.imag() == 0 );
}
inline bool operator!=(const comp& x,const comp& y){
    return ( x.real() != y.real() || x.imag() != y.imag() );
}
inline bool operator!=(const comp& x,const float& y){
    return ( x.real() != y || x.imag() != 0 );
}
inline bool operator!=(const float& x,const comp& y){
    return ( x != y.real() || y.imag() != 0 );
}

comp exp(comp z){
    return comp(exp(z.real())*cos(z.imag()),exp(z.real())*sin(z.imag()));
}

comp pow(comp z, int n){
    if(n==0) return comp(1,0);
    else if(n>0) return z*pow(z,n-1);
    return comp(0,0);
}

double abs(comp z){
    return sqrt(pow(z.real(),2)+pow(z.imag(),2));
}


comp conjugate(comp z) {
	return comp(z.re, -z.im);
}
