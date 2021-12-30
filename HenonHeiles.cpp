#include "capd/capdlib.h"
#include <vector>
#include <iostream>
#include "capd/dynset/C11Rect2Set.hpp"
#include "capd/dynsys/FirstOrderEnclosure.hpp"
#include "capd/dynsys/FirstOrderEnclosure.h"
#include "HigherOrderEnclosure.h"

using namespace capd;
using namespace dynsys;
using namespace vectalg;
using namespace std;

#include "Node.h"
#include "Var.h"

#define LOG(x) cout << (#x) << "=" << (x) << endl;

// Functions to print any vector or matrix
template<typename VectorType>
void print_vector(VectorType v){
	for(unsigned int i=0; i<v.dimension(); i++){
		cout << v[i] << " ";
	}
	cout << endl;
}

template<typename MatrixType>
void print_matrix(MatrixType m){
	for(unsigned int i=0; i<m.dimension().first; i++){
		for(unsigned int j=0; j<m.dimension().second; j++){
			cout << m[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


// Henon-Heiles class
template <typename VectorType, typename MatrixType>
class HH : public capd::dynsys::C1DynSys<MatrixType>{
public:
	Node::DAG dag; // DAG = std::vector<Node*> 
	typedef typename VectorType::ScalarType ScalarType;
	vector<vector<Var>> w1, w2, w3; 
	vector<vector<Var>> m1, m2; 
	vector<vector<Var>> coeffs_x0, coeffs, coeffs_rem;
	vector<vector<Var>> dcoeffs, coeffs_jacRem;
	
	void set_hamiltonian(interval hamiltonian){
		this->hamiltonian = hamiltonian;
	}
	
	void set_conserveEnergy(bool p){
		this->conserveEnergy = p;
	}
	
	void set_makeSymplectic(bool p){
		this->makeSymplectic = p;
	}
	
	// Functions required by capd::dynsys::C1DynSys<MatrixType>{}
	void encloseC0Map(
                const ScalarType& t,
                const VectorType& x,
                const VectorType& xx,
                VectorType& o_phi,
                VectorType& o_rem,
                VectorType& o_enc,
                MatrixType& o_jacPhi) {throw "encloseC0Map";}
		
        VectorType Remainder(const ScalarType& t, const VectorType &iv, VectorType &out_enc) {throw "Remainder";}
        VectorType enclosure(const ScalarType& t, const VectorType& x) {throw "enclosure";}
        MatrixType JacPhi(const ScalarType& t, const VectorType &iv) {throw "JacPhi"; }
        VectorType Phi(const ScalarType& t, const VectorType &iv) {throw "Phi";}
	
	
	// Constructor	
	HH(VectorType x, int order) :
        w1(order, vector<Var> (6)), w2(order+1, vector<Var> (6)), w3(order, vector<Var> (6)),
        m1(4*order, vector<Var> (6)), m2(4*(order+1), vector<Var> (6)),
        coeffs_x0(order+1, vector<Var> (4)), coeffs(order+1, vector<Var> (4)), coeffs_rem(order+2, vector<Var> (4)), 
        dcoeffs(4*(order+1), vector<Var> (4)), coeffs_jacRem(4*(order+2), vector<Var> (4)),
        order(order), x(x),
        previous_remainder_c0(4), previous_remainder_c1(4,4) {
		
		calculate_hamiltonian(this->x);	
		this->current_time = 0;
		Node::dag = &dag;
		set_zeros();
	}
	
	
	void calculate_hamiltonian(VectorType x){
		this->hamiltonian = (sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3]))/2+sqr(x[0])*x[1]-sqr(x[1])*x[1]/3;
		cout << "Hamiltonian: " << hamiltonian << endl;
	}


	void computeODECoeffs(vector<vector<Var>>& w, vector<vector<Var>>& coeffs, int order){ //x^2, y^2, 2xy, -x-2xy, y^2-y, y^2-y-x^2
                ::computeODECoeffs(w,coeffs,order);
	}


	void computeODECoeffs(vector<vector<Var>>& w, vector<vector<Var>>& coeffs, vector<vector<Var>>& M, vector<vector<Var>>& dcoeffs, int order){
		::computeODECoeffs(w,coeffs,M,dcoeffs,order);
	}


	vector<Var> taylor_poly(vector<vector<Var>>& coeffs, interval h){		
                return ::taylor_poly(coeffs,h,order);

	}

	vector<vector<Var>> taylor_poly_d(vector<vector<Var>>& dcoeffs, interval h){
                return ::taylor_poly_d(dcoeffs,h,order);
	}
	
	void conserve_energy(vector<Var>& x){
		cout << "CONSERVE ENERGY" << endl;
		Var hamiltonian = (sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3]))/Var(interval(2))+sqr(x[0])*x[1]-sqr(x[1])*x[1]/Var(interval(3));
		Node::eval(dag);
		
		cout << "Hamiltonian after step: " << hamiltonian.node->value << endl;
		hamiltonian.setValue(this->hamiltonian);
		cout << hamiltonian.node->value << endl;
		cout << hamiltonian.node->contract() << endl;	
	}
	
	
	 void make_symplectic(vector<vector<Var>>& M){
	
		int dim = M[0].size();
		vector<vector<Var>> result(dim, vector<Var> (dim));
				
		for(int i = 0; i<dim; i++){
			for(int j = 0; j<dim; j++){
				Var c = M[dim-1-dim/2][i]*M[dim-1][j];
				for(int k = 1; k<dim; k++){
					if(k<dim/2){
						Var b = M[dim-k-1-dim/2][i] * M[dim-k-1][j];
						Var d = c+b;
						c = d;
					}
					else{
						Var b = M[dim-k-1+dim/2][i] * M[dim-k-1][j];
						Var d = c-b;
						c = d;
					}	
				}
				result[i][j]=c;
			}
		}		
	   
	   Node::eval(dag);
	   
	   cout << "Original matrix" << endl;
	   for(int i=0; i<dim; i++){
		   for(int j=0; j<dim ;j++){
			   cout << M[i][j].node->value;
			   cout << " ";
		   }
		   cout << endl;
	   }
	   
	   vector<vector<interval>> values(dim, vector<interval> (dim));
	   for(int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				values[i][j] = M[i][j].node->value;
			}
	   }
	   
	   for(int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				if(j==i+dim/2){
					result[i][j].setValue(1.);
				}
				else if(i==j+dim/2){
					result[i][j].setValue(-1.);
				}
				else{
					result[i][j].setValue(0.);
				}
				cout << result[i][j].node->contract() << endl;
			}
		}
		
	   cout<< "After contraction" << endl;
	   for(int i=0; i<dim; i++){
		   for(int j=0; j<dim ;j++){
			   cout << M[i][j].node->value;
			   cout << " ";
		   }
		   cout << endl;
	   }
	   
	   
	   double diff = 0;
	   for(int i=0; i<dim; i++){
		  for(int j=0; j<dim; j++){
			  cout << width(values[i][j])-width(M[i][j].node->value)<< endl;
			  double n = (width(values[i][j])-width(M[i][j].node->value))/width(values[i][j]);
			  diff+=n;
			}
	   }
	   diff = diff/(dim*dim)*100;
	   cout << "Total percentage difference: " << diff << "%" << endl;
	}
	
	ScalarType getStep() const {return this->h;}

    void encloseC1Map(
      const ScalarType& t, 
      const VectorType& x, //x0
      const VectorType& xx, //x
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc){
		  
		  cout << "encloseC1Map" << endl;
		  cout << "Current time: " << t << endl;
		  
		  cout << "x0" << endl;
		  print_vector<IVector>(x);
		  cout << "x" << endl;
		  print_vector<IVector>(xx);
		  cout << "y" << endl;
		  print_vector<IVector>(o_phi);
		  cout << "enc" << endl;
		  print_vector<IVector>(o_enc);
		  cout << "rem" << endl;
		  print_vector<IVector>(o_rem);
		
		  cout << "jacPhi" << endl;
		  print_matrix<IMatrix>(o_jacPhi);
		  cout << "jacEnc" << endl;
		  print_matrix<IMatrix>(o_jacEnc);
		  cout << "jacRem" << endl;
		  print_matrix<IMatrix>(o_jacRem);
			  
		  set_zeros();
			  
		  for(int i=0; i<4; i++){
				coeffs[0][i]=Var(xx[i]);
		  }
				
		  for(int i=0; i<4; i++){
				dcoeffs[i][i] = Var(interval(1));
		  }
		  	
		  	
		  // HOE	
		  o_enc = this->C0_enclose_1(xx, this->current_time, this->h); // enclosure for C0 part
		  cout << o_enc << endl;
		  o_jacEnc = this->C1_enclose_1(xx, this->current_time, this->h, o_enc); // enclosure for C1 part
		  cout << o_jacEnc << endl;
		  
		  
		  // FOE
		/* o_enc = this->C0_enclose(xx, this->current_time, this->h); // enclosure for C0 part
		  cout << o_enc << endl;
		  o_jacEnc = this->C1_enclose(xx, this->current_time, this->h); 
		  cout << o_jacEnc << endl;*/
		  
		  
		  vector<vector<Var>> o_rem_var(1, vector<Var> (4)); // Var of o_rem
		  vector<vector<Var>> o_jacRem_var(4, vector<Var> (4)); // Var of o_jacRem
			
		  // o_rem and o_jacRem
		  for(int i=0; i<4; i++){
				coeffs_rem[0][i] = Var(o_enc[i]);
		  }
		  
		  for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				coeffs_jacRem[i][j]=Var(o_jacEnc[i][j]);
			}
		  }
	
		this->computeODECoeffs(this->w2, this->coeffs_rem, this->m2, this->coeffs_jacRem, this->order+1); 
		
		// h^(order+1)
		double h_ = 1;
		for(int i=0; i<order+1; i++){
			h_ *= h;
		}
		
		Var hh = Var(interval(h_));
		
		for(int i=0; i<4; i++){
			o_rem_var[0][i] = this->coeffs_rem[order+1][i]*hh;
		}
		
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				o_jacRem_var[i][j] = this->coeffs_jacRem[4*(order+1)+i][j]*hh;
			}
		}
		
		// calculate o_phi and o_jacPhi
		this->computeODECoeffs(this->w1, this->coeffs, this->m1, this->dcoeffs, this->order); 
		vector<vector<Var>> o_jacPhi_var = taylor_poly_d(this->dcoeffs, this->h);

		// calculate o_phi
		for(int i=0; i<4; i++){
			this->coeffs_x0[0][i]=Var(x[i]);
		}
		
		this->computeODECoeffs(this->w3, this->coeffs_x0, this->order);
		vector<Var> o_phi_var = taylor_poly(this->coeffs_x0, this->h);
		
		
	
	
		
		//**************************************************************
		//conserve energy
		vector<Var> tmp(4);
		for(int i=0; i<4; i++){
			tmp[i] = Var(interval(0));
		}
		
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				tmp[i] += o_jacPhi_var[i][j]*(coeffs[0][j]-coeffs_x0[0][j]); //jacPhi*(x-x0)
			}
		}
		
		vector<Var> flow(4);
		for(int i=0; i<4; i++){
			flow[i] = o_phi_var[i]+tmp[i]+o_rem_var[0][i];
		}	
		
		Var hamiltonian = (sqr(flow[0])+sqr(flow[1])+sqr(flow[2])+sqr(flow[3]))/Var(interval(2))+sqr(flow[0])*flow[1]-sqr(flow[1])*flow[1]/Var(interval(3));
	
		
		
		
		 //**************************************************************
		// make symplectic
		vector<vector<Var>> sympl(4, vector<Var> (4));
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				sympl[i][j] = o_jacPhi_var[i][j]+o_jacRem_var[i][j];
			}
		}
		
		int dim = 4;
		vector<vector<Var>> result(dim, vector<Var> (dim));
				
		for(int i = 0; i<dim; i++){
			for(int j = 0; j<dim; j++){
				Var c = sympl[dim-1-dim/2][i]*sympl[dim-1][j];
				for(int k = 1; k<dim; k++){
					if(k<dim/2){
						Var b = sympl[dim-k-1-dim/2][i] * sympl[dim-k-1][j];
						Var d = c+b;
						c = d;
					}
					else{
						Var b = sympl[dim-k-1+dim/2][i] * sympl[dim-k-1][j];
						Var d = c-b;
						c = d;
					}	
				}
				result[i][j]=c;
			}
		}		
	   
	   	Node::eval(dag);
		
		cout << "Hamiltonian after step: " << hamiltonian.node->value << endl;
		hamiltonian.setValue(this->hamiltonian);
		cout << hamiltonian.node->value << endl;
		
		if(this->conserveEnergy){
			
		cout << hamiltonian.node->contract() << endl;
		} 
	   
	
		
		if(this->makeSymplectic){
			  for(int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				if(j==i+dim/2){
					result[i][j].setValue(1.);
				}
				else if(i==j+dim/2){
					result[i][j].setValue(-1.);
				}
				else{
					result[i][j].setValue(0.);
				}
				cout << result[i][j].node->contract() << endl;
			}
		}
		}
		
	

		// copy values
		for(int i=0; i<4; i++){
			o_rem[i] = o_rem_var[0][i].node->value;
			o_phi[i] = o_phi_var[i].node->value;
		}
				
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				o_jacRem[i][j] = o_jacRem_var[i][j].node->value;
				o_jacPhi[i][j] = o_jacPhi_var[i][j].node->value;
			}
		}
		
		Node::dag->clear();
		this->current_time += this->h;
	}
	
	
	VectorType C0_enclose_1(VectorType u, double current_time, double time_step){
		VectorType enc = high_C0_enclose_1(u, current_time, time_step, this->order, this->previous_remainder_c0);
		return enc;
	}

	MatrixType C1_enclose_1(VectorType u, double current_time, double time_step, VectorType enclosure_c0){
		cout << "Liczymy enclosure C1" << endl;
		MatrixType enc = high_C1_enclose_1(u, current_time, time_step, this->order, enclosure_c0, this->previous_remainder_c1);
		return enc;
	}

	static VectorType C0_enclose(VectorType u, double current_time, double time_step){
		VectorType enc = FirstOrderEnclosure::enclosure(HenonHeilesVectorField,
													    current_time,
													    u,
													    time_step);										    
		//cout << enc << endl;
		return enc;
	}
	
	
	
	static MatrixType C1_enclose(VectorType u, double current_time, double time_step){
		VectorType enc = C0_enclose(u, current_time, time_step);
		MatrixType jacEnc = FirstOrderEnclosure::jacEnclosure(HenonHeilesVectorField,
														      current_time,
														      time_step,
														      enc,
														      EuclLNorm<VectorType, MatrixType>());
		//cout << jacEnc << endl;
		return jacEnc;
	}
	
	void set_timestep(double timestep){
		this->h=timestep;
	}

private:
        double h;
        interval hamiltonian;
        bool conserveEnergy;
        bool makeSymplectic;
        double current_time;
        static IMap HenonHeilesVectorField;

        int order;
        VectorType x;
        VectorType previous_remainder_c0;
        MatrixType previous_remainder_c1;

        // Fills w1, w2, w3, m1, m2, coeffs_x0, coeffs, coeffs_rem, dcoeffs, coeffs_jacRem with zeros
        void set_zeros();
        void set_zeros(vector<Var>& v){	static Var zero= Var(interval(0.)); for(auto& u: v) u = zero; }
        void set_zeros(vector<vector<Var>>& v){	for(auto& u: v) set_zeros(u); }
}; // class HH

template <typename VectorType, typename MatrixType>
IMap HH<VectorType, MatrixType>::HenonHeilesVectorField("var:x,y,px,py;fun:px,py,-x-2*x*y,y^2-y-x^2;");

template <typename VectorType, typename MatrixType>
void HH<VectorType, MatrixType>::set_zeros(){

        // Fill w1, w2, w3, m1, m2 with zeros
        set_zeros(w1);
        set_zeros(w2);
        set_zeros(w3);
        set_zeros(m1);
        set_zeros(m2);
        
        // Fill coeffs_x0, coeffs, coeffs_rem, dcoeffs, coeffs_jacRem with zeros
        set_zeros(coeffs);
        set_zeros(coeffs_x0);
        set_zeros(coeffs_rem);
        set_zeros(coeffs_rem);
        
        set_zeros(dcoeffs);
        set_zeros(coeffs_jacRem);
}



void test_C0(){

	 int order = 3;
	 IVector u = IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1; //x,y,px,py
     //IVector u = IVector({1,2,3,4}); //x,y,px,py

     HH <IVector, IMatrix> hh(u, order);
     for(int i=0; i<4; i++){
		 hh.coeffs[0][i] = Var(u[i]);
	 }
     hh.computeODECoeffs(hh.w1, hh.coeffs, order);
     Node::eval(hh.dag);
     
     cout << "Solution" << endl;
     for(int i=0; i<order+1; i++){
		 for(int j=0; j<4; j++){
			 cout << hh.coeffs[i][j].node->value << " ";
		 }
		 cout << endl;
	 }

     // IMap
	 IVector* c = new IVector[order+1];
	 c[0] = u;
	 IVector zeros = IVector({0,0,0,0});
	 for(int i=1; i<order+1; i++){
		c[i]=zeros;
	 }

	 IMap HenonHeilesVectorField("var:x,y,px,py;fun:px,py,-x-2*x*y,y^2-y-x^2;");
	 IOdeSolver solver(HenonHeilesVectorField, 20);
	 HenonHeilesVectorField.computeODECoefficients(c, order);

	 cout << "Solution IMap" << endl;
	 for(int i=0; i<order+1; i++){
		  cout << c[i] << endl;
	 }
}

void test_C1(){

	int order = 3;
	IVector u = IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1; //x,y,px,py
  //IVector u = IVector({1,2,3,4}); //x,y,px,py

    HH <IVector, IMatrix> hh(u, order);
    for(int i=0; i<4; i++){
		 hh.coeffs[0][i] = Var(u[i]);
		 hh.dcoeffs[i][i] = Var(interval(1.));
	 }
	 
    hh.computeODECoeffs(hh.w1, hh.coeffs, hh.m1, hh.dcoeffs, order);
    Node::eval(hh.dag);
    
    cout << "Solution" << endl;
     for(int i=0; i<4*(order+1); i++){
			 for(int j=0; j<4; j++){
				 cout << hh.dcoeffs[i][j].node->value << " ";
			 }
			 cout << endl;
	 }
    

	// IMap
	IVector zeros = IVector({0,0,0,0});
	IVector* m = new IVector[order+1];
	m[0]=u;
	for(int i=1; i<order+1; i++){
		m[i]=zeros;
	}

	IMatrix* M = new IMatrix[order+1];
	IMatrix Q(4,4);
	for(int i=0; i<order+1; i++){
		M[i]=Q;
	}

	Q(1,1) = Q(2,2) = Q(3,3) = Q(4,4) = 1;
	M[0]=Q;
	IMap HenonHeilesVectorField("var:x,y,px,py;fun:px,py,-x-2*x*y,y^2-y-x^2;");
	IOdeSolver solver(HenonHeilesVectorField, 20);
	HenonHeilesVectorField.computeODECoefficients(m,M,order);
	cout << "IMap solution" << endl;
	for(int i=0; i<order+1; i++){
		cout << M[i] << endl;
	}
}

void test_enclosure(){

	IVector x = IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1;
	double current_time = 0;
	double time_step = 0.1;

	IVector enc = HH<IVector, IMatrix>::C0_enclose(x, current_time, time_step);
	cout << enc << endl;

	IMatrix jacEnc = HH<IVector, IMatrix>::C1_enclose(x, current_time, time_step);
	cout << jacEnc << endl;
}



void test_encloseC1Map(){
	int dim = 4;
	IVector y(dim), rem(dim), enc(dim);
	IMatrix jacPhi(dim,dim), jacEnc(dim,dim), jacRem(dim,dim);
	
	int order = 10;
	double h = 0.125;
	IVector x = IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1;
	IVector x0 = IVector({1,2,3,4});

	HH <IVector, IMatrix> henon(x, order);
	henon.set_timestep(h);
	henon.set_makeSymplectic(false);
	henon.set_conserveEnergy(false);
	henon.encloseC1Map(0,x0,x,y,rem,enc,jacPhi,jacRem,jacEnc);
	
	cout << "TEST" << endl;
	
	cout << "y" << endl;
	print_vector<IVector>(y);
	cout << "enc" << endl;
	print_vector<IVector>(enc);
	cout << "rem" << endl;
	print_vector<IVector>(rem);
	
	cout << "jacPhi" << endl;
	print_matrix<IMatrix>(jacPhi);
	cout << "jacEnc" << endl;
	print_matrix<IMatrix>(jacEnc);
	cout << "jacRem" << endl;
	print_matrix<IMatrix>(jacRem);
	
	vector<interval> tmp(4);
	for(int i=0; i<4; i++){
		tmp[i] = interval(0);
	}
		
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			tmp[i] += jacPhi[i][j]*(x[j]-x0[j]);
		}
	}
		
	vector<interval> flow(4);
	for(int i=0; i<4; i++){
		flow[i] = y[i]+tmp[i]+rem[0];
		cout << flow[i] << " ";
	}
}

void test_symplectic(){
	
	int order = 10;
	double h = 0.0136719;
	IVector x = IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1;
	IVector x0 = IVector({1,2,3,4});

	HH <IVector, IMatrix> henon(x, order);
	vector<vector<Var>> s(2, vector<Var> (2)); 
	s[0][0]=Var(interval(0,2));
	s[0][1]=Var(interval(-1, 1));
	s[1][0]=Var(interval(-1,-0.5));
	s[1][1]=Var(interval(-3,-2));
	
	//henon.make_symplectic(s);
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			cout << s[i][j].node->value << " ";
		}
		cout << endl;
	}
}

void test_1(bool symplectic, bool hamiltonian){
	/*
	 * Ustalamy x=0, py=0, y oraz H, wyliczamy px
	 */
	
	vector<vector<interval>> set1_res(100, vector<interval>(4));
	vector<vector<interval>> set2_res(100, vector<interval>(4));
	
	vector<vector<interval>> set1_res_jacPhi(400, vector<interval>(4));
	vector<vector<interval>> set2_res_jacPhi(400, vector<interval>(4));
	
	int n = 100;
	int order = 5;
	double h = 0.125;
	
	interval H(1./6);
	interval y(0.3407722818724012);
	interval py(0);
	interval x(0);
	
	interval px = sqrt(2*H-sqr(py)+sqr(y)*(y*2./3-1));
	IVector v = IVector({x,y,px,py});

	HH <IVector, IMatrix> henon1(v, order);
	henon1.set_timestep(h);
	henon1.set_hamiltonian(H);
	henon1.set_makeSymplectic(false);
	henon1.set_conserveEnergy(false);
	
	C1Rect2Set set1(v);
	for(int i=0; i<n; i++){
		set1.move(henon1);
		for(int j=0; j<4; j++){
			set1_res[i][j] = IVector(set1)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set1_res_jacPhi[4*i+k][j] = IMatrix(set1)[k][j];
			}
		}
	}
	
	HH <IVector, IMatrix> henon2(v, order);
	henon2.set_timestep(h);
	henon2.set_hamiltonian(H);
	henon2.set_makeSymplectic(symplectic);
	henon2.set_conserveEnergy(hamiltonian);
	
	C1Rect2Set set2(v);
	for(int i=0; i<n; i++){
		set2.move(henon2);
		for(int j=0; j<4; j++){
			set2_res[i][j] = IVector(set2)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set2_res_jacPhi[4*i+k][j] = IMatrix(set2)[k][j];
			}
		}
	}
	
	for(int i=0; i<n; i++){
		double diff = 0;
		double diff_jacPhi = 0;
		for(int j=0; j<4; j++){
			double n = (width(set1_res[i][j])-width(set2_res[i][j]))/width(set1_res[i][j]);
			diff+=n;
		}
		cout << "Iteration: " << i << " difference: " << diff/4.*100 << "%" << endl;
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				double n = (width(set1_res_jacPhi[4*i+k][j])-width(set2_res_jacPhi[4*i+k][j]))/width(set1_res_jacPhi[4*i+k][j]);
				diff_jacPhi+=n;
			}
		}
		cout << "Iteration: " << i << " difference: " << diff_jacPhi/16.*100 << "%" << endl;
	}
	
	LOG(IVector(set1));
	LOG(IMatrix(set1));
	LOG(IVector(set2));
	LOG(IMatrix(set2));
}

//Podajemy tylko order, h, x, czy sympl czy hamiltonian

void test_2(int order, double h, int n, IVector x, bool symplectic, bool hamiltonian){
	
	vector<vector<interval>> set1_res(n, vector<interval>(4));
	vector<vector<interval>> set2_res(n, vector<interval>(4));
	
	vector<vector<interval>> set1_res_jacPhi(4*n, vector<interval>(4));
	vector<vector<interval>> set2_res_jacPhi(4*n, vector<interval>(4));

	HH <IVector, IMatrix> henon1(x, order);
	henon1.set_timestep(h);
	henon1.set_makeSymplectic(false);
	henon1.set_conserveEnergy(false);
	
	C1Rect2Set set1(x);
	for(int i=0; i<n; i++){
		set1.move(henon1);
		LOG(IVector(set1));
		LOG(IMatrix(set1));
		for(int j=0; j<4; j++){
			set1_res[i][j] = IVector(set1)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set1_res_jacPhi[4*i+k][j] = IMatrix(set1)[k][j];
			}
		}
	}
	
	HH <IVector, IMatrix> henon2(x, order);
	henon2.set_timestep(h);
	henon2.set_makeSymplectic(symplectic);
	henon2.set_conserveEnergy(hamiltonian);
	
	C1Rect2Set set2(x);
	for(int i=0; i<n; i++){
		set2.move(henon2);
		for(int j=0; j<4; j++){
			set2_res[i][j] = IVector(set2)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set2_res_jacPhi[4*i+k][j] = IMatrix(set2)[k][j];
			}
		}
	}
	
	for(int i=0; i<n; i++){
		double diff = 0;
		double diff_jacPhi = 0;
		for(int j=0; j<4; j++){
			double n = (width(set1_res[i][j])-width(set2_res[i][j]))/width(set1_res[i][j]);
			diff+=n;
		}
		cout << "Iteration: " << i << " difference: " << diff/4.*100 << "%" << endl;
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				double n = (width(set1_res_jacPhi[4*i+k][j])-width(set2_res_jacPhi[4*i+k][j]))/width(set1_res_jacPhi[4*i+k][j]);
				diff_jacPhi+=n;
			}
		}
		cout << "Iteration: " << i << " difference: " << diff_jacPhi/16.*100 << "%" << endl;
	}
	
	LOG(IVector(set1));
	LOG(IMatrix(set1));
	LOG(IVector(set2));
	LOG(IMatrix(set2));
}



void test_3(int order, double h, int n, bool symplectic, bool hamiltonian){
	vector<vector<interval>> set1_res(n, vector<interval>(4));
	vector<vector<interval>> set2_res(n, vector<interval>(4));
	
	vector<vector<interval>> set1_res_jacPhi(4*n, vector<interval>(4));
	vector<vector<interval>> set2_res_jacPhi(4*n, vector<interval>(4));
	
	interval y(0.3407722818724012);
	interval py(0);
	interval x(0);
	interval H(1./6);
	
	interval px = sqrt(2*H-sqr(y)-sqr(py)+sqr(y)*y*2./3);
	IVector x0({x,y,px,py});
	
	interval delta = interval(-1,1)*1e-2;
	y += delta;
	py += delta;
	px = sqrt(2*H-sqr(y)-sqr(py)+sqr(y)*y*2./3);
	
	IMatrix C(4,4);
	C(2,2) = C(4,4) = interval(1);
	C(3,2) = (-1)*(y+sqr(y))/px;
	C(3,4) = (-1)*py/px;
	
	IVector r({0,delta,0,delta});
	C1Rect2Set set1(x0,C,r);
	C1Rect2Set set2(x0,C,r);

	HH <IVector, IMatrix> henon1(x0, order);
	henon1.set_timestep(h);
	henon1.set_hamiltonian(H);
	henon1.set_makeSymplectic(false);
	henon1.set_conserveEnergy(false);
	
	for(int i=0; i<n; i++){
		set1.move(henon1);
		for(int j=0; j<4; j++){
			set1_res[i][j] = IVector(set1)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set1_res_jacPhi[4*i+k][j] = IMatrix(set1)[k][j];
			}
		}
	}
	
	HH <IVector, IMatrix> henon2(x0, order);
	henon2.set_timestep(h);
	henon2.set_hamiltonian(H);
	henon2.set_makeSymplectic(symplectic);
	henon2.set_conserveEnergy(hamiltonian);
	
	for(int i=0; i<n; i++){
		set2.move(henon2);
		for(int j=0; j<4; j++){
			set2_res[i][j] = IVector(set2)[j];
		}
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				set2_res_jacPhi[4*i+k][j] = IMatrix(set2)[k][j];
			}
		}
	}
	
	for(int i=0; i<n; i++){
		double diff = 0;
		double diff_jacPhi = 0;
		for(int j=0; j<4; j++){
			double n = (width(set1_res[i][j])-width(set2_res[i][j]))/width(set1_res[i][j]);
			diff+=n;
		}
		cout << "Iteration: " << i << " difference: " << diff/4.*100 << "%" << endl;
		
		for(int k=0; k<4; k++){
			for(int j=0; j<4; j++){
				double n = (width(set1_res_jacPhi[4*i+k][j])-width(set2_res_jacPhi[4*i+k][j]))/width(set1_res_jacPhi[4*i+k][j]);
				diff_jacPhi+=n;
			}
		}
		cout << "Iteration: " << i << " difference: " << diff_jacPhi/16.*100 << "%" << endl;
	}
	
	LOG(IVector(set1));
	LOG(IMatrix(set1));
	LOG(IVector(set2));
	LOG(IMatrix(set2));
}


int main () {
	cout.precision(10);
        try{
                //test_enclosure();
                //test_C0();
                //test_C1();
                //test_symplectic();
                //test_encloseC1Map();
                
                //test_1(true, true);
                
                //test_2(5, 0.013, 100, IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1, false, true);
                //test_2(5, 0.05, 29, IVector({1,2,3,4}), true, true);
                //test_2(5, 0.05, 29, IVector({interval(0.5,1.5),2,3,4}), true, true); 1%
                //test_2(5, 0.125, 100, IVector({0.1,0.2,0,0}), true, true); 
                //test_2(5, 0.25, 40, IVector({0.,0.2,0.,0.})+IVector({1,1,1,1})*interval(-1,1)*1e-1, true, true); // tu  2,7% i 4,2%
                test_2(5, 0.25, 13, IVector({0.,0.2,0.,0.})+IVector({1,1,1,1})*interval(-1,1)*1e-1, true, true); // tu 6,08, 10,66
                //test_2(5, 0.05, 80, IVector({0,0.5,0,0}), true, true);
                //test_2(10, 1./32, 32, IVector({1,2,3,4})+IVector({1,1,1,1})*interval(-1,1)*1e-1, true, true);
                
                //test_3(5, 1./16, 60, true, true); // 0.02%
        }catch(std::exception& e){
                cout << "Exception: " << e.what() << endl;
        }catch(const char* e){
                cout << "Exception: " << e << endl;
        }
	return 0;
}
