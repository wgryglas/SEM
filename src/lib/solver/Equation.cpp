
//#include"Equation.h"


//namespace SEM{ namespace las
//{
//    Equation::Equation(ScalarField& field)
//        :field(field), sparseBuild(false)
//    {
//        initialSol.setZero(field.Tacc.size());
//        for(int i=0;i<field.Tacc.size();i++)
//            initialSol[i]=field.Tacc[i];
//    }

//    void Equation::applyDirichletBoundaryCondition()
//    {
//        //Apply Dirichlet BC--> set at diag. of coresponding equation 1, and in coresponding row in rhs value of dirichlet BC
//        for(int row=0; row<field.dirichletMask.size();row++)
//        {
//            if(field.dirichletMask[row])
//            {
//                coeffs.push_back( MCoeff(row,row,1) );
//                rhsVector[row]= field.Tacc[row];
//            }
//        }
//    }

//    void Equation::buildSparseMatrix()
//    {
//        matrix.resize(field.Tacc.size(),field.Tacc.size());
//        matrix.setFromTriplets(coeffs.begin(),coeffs.end());
//        sparseBuild = true;
//    }

//    inline bool Equation::isBuild()
//    {
//        return sparseBuild;
//    }

//    inline void Equation::setAsUnBuild()
//    {
//        this->sparseBuild = false;
//    }

//    Equation& Equation::operator+(Equation& eq)
//    {
//        if(isBuild() && eq.isBuild())
//        {
//            this->matrix+=eq.matrix;
//        }
//        else
//        {
//            coeffs.reserve(coeffs.size()+eq.coeffs.size());
//            for(vector<MCoeff>::iterator i = eq.coeffs.begin();i!=eq.coeffs.end();i++)
//                coeffs.push_back(*i);

//            setAsUnBuild();
//        }

//        this->rhsVector+=eq.rhsVector;

//        return *this;
//    }

//    Equation& Equation::operator-(Equation& eq)
//    {
//        if(isBuild() && eq.isBuild())
//        {
//            this->matrix-=eq.matrix;
//        }
//        else
//        {
//            coeffs.reserve(coeffs.size()+eq.coeffs.size());
//            for(vector<MCoeff>::iterator i = eq.coeffs.begin();i!=eq.coeffs.end();i++)
//                coeffs.push_back(MCoeff(i->row(),i->col(),-(i->value())));

//            setAsUnBuild();
//        }

//        this->rhsVector-=eq.rhsVector;
//        return *this;
//    }

//    Equation& Equation::operator=(Vector& v)
//    {
//        rhsVector = rhsVector + v;
//        return *this;
//    }

//    ostream& operator<<(ostream& o,Equation& e)
//    {
//        cout<<"Matirx::"<<endl;
//        for(int i=0;i<e.matrix.rows();i++)
//        {
//            for(int j=0;j<e.matrix.cols();j++)
//                cout<<e.matrix.coeff(i,j)<<"|";
//            cout<<endl;
//        }

//        cout<<"Vector::"<<endl;
//        for(int i=0; i<e.rhsVector.size(); i++)
//            cout<<e.rhsVector[i]<<endl;

//        return o;
//    }
//}//las
//}//SEM
