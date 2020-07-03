#ifndef _Position_
#define _Position_

class Position {

private:

protected:
    double m_x, m_y, m_z;

public:
    // constructors
    Position();
    Position(double x,double y, double z);
    // destructor
    ~Position();
    // methods
    void SetX(double);
    void SetY(double);
    void SetZ(double);
    double GetX();
    double GetY();
    double GetZ();
    double GetDist();
};

#endif 
