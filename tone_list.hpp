#pragma once

class Tone{
public:
    double t[88];
    
    Tone(){
        for(int i=0;i<88;i++){
            t[i]=440*pow(2,(i-48)/12.0);
        }
    }
    
    int searchTone(double f){
        return search(f,0,87);
    }
    
    int search(double f, int min, int max){
        if(min>max || f<t[min] || f>t[max]){
            printf("No corresponging element\n");
            return -1;
        }
        
        if(min+1==max){
            if(f-t[min]>t[max]-f) return max;
            else return min;
        }
        else{
            int i_mid=(max+min)/2;
            if(t[i_mid]>f){
                return search(f,min,i_mid);
            }
            else{
                return search(f,i_mid,max);
            }
        }
    }
    
    string toneName(int i){
        string No, Al;
        if(0<=i && i<=2){
            No="0";
            if(i==0) Al="A";
            if(i==1) Al="A#";
            if(i==2) Al="B";
            string ret=Al+No;
            return ret;
        }
        else if(3<=i && i<=14){
            No="1";
            i-=3;
        }
        else if(15<=i && i<=26){
            No="2";
            i-=15;
        }
        else if(27<=i && i<=38){
            No="3";
            i-=27;
        }
        else if(39<=i && i<=50){
            No="4";
            i-=39;
        }
        else if(51<=i && i<=62){
            No="5";
            i-=51;
        }
        else if(63<=i && i<=74){
            No="6";
            i-=63;
        }
        else if(75<=i && i<=86){
            No="7";
            i-=75;
        }
        else if(i==87){
            No="8";
            Al="C";
            string ret=Al+No;
            return ret;
        }
        
        switch(i){
            case 0:
                Al="C";
                break;
            case 1:
                Al="C#";
                break;
            case 2:
                Al="D";
                break;
            case 3:
                Al="D#";
                break;
            case 4:
                Al="E";
                break;
            case 5:
                Al="F";
                break;
            case 6:
                Al="F#";
                break;
            case 7:
                Al="G#";
                break;
            case 8:
                Al="A";
                break;
            case 9:
                Al="A#";
                break;
            case 10:
                Al="B";
                break;
            case 11:
                Al="B#";
                break;
        }
        string ret=Al+No;
        return ret;
    }
};
