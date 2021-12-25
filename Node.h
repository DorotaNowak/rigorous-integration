struct Node{
  enum Result {EMPTY=-1,NOIMPROVEMENT=0,CONTRACTION=1};
  explicit Node(capd::interval p) 
    : value(p), left(nullptr), right(nullptr)
  {}
  
  Node(Node* left, Node* right)
    : left(left), right(right)
  {}
  
  virtual void eval() {}
  virtual Result contract() { return NOIMPROVEMENT; }
  inline void setValue(capd::interval p) { this->value = p; this->isnan = false; }

  Node* left;
  Node* right;
  capd::interval value;
  bool isnan = false;
  bool isOutNode = false;
  
  typedef std::vector<Node*> DAG;
  static DAG* dag;  
  static void eval(DAG& dag){
    for(auto* n : dag) {
      n->eval();
      if(n->isOutNode and !n->isnan and !isSingular(n->value)) return;
    }
  }
};

Node::DAG* Node::dag = nullptr;

struct VarNode : public Node{
  VarNode() : Node(nullptr,nullptr) {}
  void eval(){ }
  virtual Result contract() { return NOIMPROVEMENT; }
};

struct ConstNode : public Node{
  explicit ConstNode(capd::interval c) : Node(c), c(c) {}
  void eval(){ this->value = this->c; }

  const capd::interval c;
};

template<class BasicOp>
struct BinaryNode : public Node, BasicOp{
  using Node::Node;
  void eval(){ 
    if(!left->isnan and !right->isnan)
      BasicOp::eval(this,left, right);
    else 
      this->isnan = true;
  }

  Result contract(){
    Result c = NOIMPROVEMENT;
    while(true){
      capd::interval L = BasicOp::getL(left->value,right->value,this->value);
      capd::interval R = BasicOp::getR(left->value,right->value,this->value);
      bool r = L.leftBound()  > left->value.leftBound() or
               L.rightBound() < left->value.rightBound() or
               R.leftBound()  > right->value.leftBound() or
               R.rightBound() < right->value.rightBound();
      if(!intersection(L,left->value,left->value) or !intersection(R,right->value,right->value)) return EMPTY;
      if(r) c=CONTRACTION;
      else break;
    }
    Result l,r; 
    if( (l=left->contract())==EMPTY or (r=right->contract())==EMPTY) return EMPTY;
    return Result(c or l or r);
  }
};

template<class BasicOp>
struct NodeWithExplicitBound : public BasicOp{
  inline void eval(Node* z, Node* l, Node* r) { 
    BasicOp::eval(z,l,r);
    if(!intersection(z->value,bound,z->value))
      z->isnan = true;
  } 
  using BasicOp::getL;
  using BasicOp::getR;
  capd::interval bound;
};

struct AddOp{
  inline static void eval(Node* z, Node* l, Node* r) { z->value =  l->value + r->value; z->isnan = false; } 
  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { return z-r;} 
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { return z-l;} 
};

struct SubOp{
  inline static void eval(Node* z, Node* l, Node* r) { z->value =  l->value - r->value; z->isnan = false; } 
  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { return z+r;} 
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { return l-z;} 
};

struct MulOp{
  inline static void eval(Node* z, Node* l, Node* r) { z->value =  l->value * r->value; z->isnan = false; } 
  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { return !isSingular(r) ? z/r : l; }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { return !isSingular(l) ? z/l : r; } 
};

struct DivOp{
  inline static void eval(Node* z, Node* l, Node* r) { 
    if(!isSingular(r->value)){
      z->value = l->value / r->value; 
      z->isnan = false;
    } else
      z->isnan = true;
  } 
  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { return z*r; }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { return !isSingular(z) ? l/z : r; } 
};

struct XByNorm3Op{

  inline static interval f(interval x, interval y) { 
    return x/power(sqrt(sqr(x)+sqr(y)),3);
  }  
  inline static interval getLeftIntersectionPoint(interval p){
    const static interval sqrt2 = sqrt(interval(2.0));
    return sqrt2*p.leftBound();
  }
  inline static interval getRightIntersectionPoint(interval p){
    const static interval sqrt2 = sqrt(interval(2.0));
    return sqrt2*p.rightBound();
  }

  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { 
    if(!intersection(l,z*power(sqrt(sqr(l)+sqr(r)),3),z)) return l.rightBound()+1;
    return z;
  }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { 
    if(isSingular(z))
      return r;
    // z = l/(l^2+r^2)^(3/2)
    // l/z = (l^2+r^2)^(3/2)
    // (l/z)^(2/3) = l^2+r^2
    // r^2 = (l/z)^(2/3) - l^2
    capd::interval u = l/z;
    u = power(u,interval(2)/3)-sqr(l);
    if(u.rightBound()<0) return r.rightBound()+1; // equation is inconsistent, so return something which has empty intersection with r
    if(u.leftBound()<0) u.setLeftBound(0.);
    if(r>0){
      intersection(r,sqrt(u),r);
    }else if(r<0){
      intersection(r,-sqrt(u),r);
    } else{
      intersection(r,interval(-1,1)*sqrt(u),r);
    }
    return r;
  } 
};

struct XByNorm2Op{

  inline static interval f(interval x, interval y) { 
    return x/(sqr(x)+sqr(y));
  }  
  inline static interval getLeftIntersectionPoint(interval p){
    return p.leftBound();
  }
  inline static interval getRightIntersectionPoint(interval p){
    return p.rightBound();
  }

  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { 
    if(!intersection(l,z*(sqr(l)+sqr(r)),z)) return l.rightBound()+1;
    return z;
  }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { 
    if(isSingular(z))
      return r;
    // z = l/(l^2+r^2)
    // l/z = (l^2+r^2)
    // r^2 = l/z - l^2= l(1/z-l)
    capd::interval u;
    intersection(l*(1./z - l),l/z-sqr(l),u);
    if(u.rightBound()<0) return r.rightBound()+1; // equation is inconsistent, so return something which has empty intersection with r
    if(u.leftBound()<0) u.setLeftBound(0.);
    if(r>0){
      intersection(r,sqrt(u),r);
    }else if(r<0){
      intersection(r,-sqrt(u),r);
    } else{
      intersection(r,interval(-1,1)*sqrt(u),r);
    }
    return r;
  } 
};

struct X2ByNorm5Op{

  inline static interval f(interval x, interval y) { 
    interval z = sqr(x);
    return z/power(sqrt(z+sqr(y)),5);
  }  
  inline static interval getLeftIntersectionPoint(interval p){
    const static interval sqrtOneHalf = sqrt(interval(1.5));
    return sqrtOneHalf*p.leftBound();
  }
  inline static interval getRightIntersectionPoint(interval p){
    const static interval sqrtOneHalf = sqrt(interval(1.5));
    return sqrtOneHalf*p.rightBound();
  }

  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { 
    return l;
    if(z.rightBound()<0) return l.rightBound()+1; //whatever which has empty intersection with l
    //  l^2 = (z^2*
    interval u = sqrt(power(sqr(z)*(sqr(l)+sqr(r)),5));
    if(l>0){
      if( !intersection(l,u,z) ) return l.rightBound()+1;
      return z;
    }
    if(l<0){
      if( !intersection(l,-u,z) ) return l.rightBound()+1;
      return z;    
    }
    if( !intersection(l,interval(-1,1)*u,z) ) return l.rightBound()+1;
    return z;    
  }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { 
    if(isSingular(z))
      return r;
    // z = l^2/(l^2+r^2)^(5/2)
    // l^2/z = (l^2+r^2)^(5/2)
    // (l^2/z)^(2/5) = l^2+r^2
    // r^2 = (l^2/z)^(2/5) - l^2
    interval l2 = sqr(l);
    capd::interval u = l2/z;
    u = power(u,interval(2)/5)-l2;
    if(u.rightBound()<0) return r.rightBound()+1; // equation is inconsistent, so return something which has empty intersection with r
    if(u.leftBound()<0) u.setLeftBound(0.);
    if(r>0){
      intersection(r,sqrt(u),r);
    }else if(r<0){
      intersection(r,-sqrt(u),r);
    } else{
      intersection(r,interval(-1,1)*sqrt(u),r);
    }
    return r;
  } 
};

/**
 * It represents expression of the form
 *    L^q/sqrt(L^2+R^2)^b
 */ 
template<class BasicOp>
struct XByNormOp : public BasicOp{
  
  inline static bool ovlap(interval x, interval y){
    return x.rightBound()>=y.leftBound() and y.rightBound()>=x.leftBound();
  }
  inline static void eval(Node* z, Node* l, Node* r) { 
    const interval& p = l->value;
    const interval& q = r->value;
    
    bool L = isSingular(p), R = isSingular(q);
    if(L and R) {
      z->isnan = true;
      return;
    }
  
    // take corners
    z->value = BasicOp::f(p.leftBound(),q.leftBound());
    z->value = intervalHull(z->value,BasicOp::f(p.leftBound(),q.rightBound()));
    z->value = intervalHull(z->value,BasicOp::f(p.rightBound(),q.rightBound()));
    z->value = intervalHull(z->value,BasicOp::f(p.rightBound(),q.leftBound()));
    
    // check for zeroes in x and y
    if(L){
      z->value = intervalHull(z->value,BasicOp::f(0.,q.leftBound()));
      z->value = intervalHull(z->value,BasicOp::f(0.,q.rightBound()));
    }
    if(R){
      z->value = intervalHull(z->value,BasicOp::f(p.leftBound(),0.));
      z->value = intervalHull(z->value,BasicOp::f(p.rightBound(),0.));
    }
    
    // check possible eight intersection points of lines 
    interval y = BasicOp::getLeftIntersectionPoint(p);
    if(ovlap(y,q))
      z->value = intervalHull(z->value,BasicOp::f(p.leftBound(),y));
    y = -y; 
    if(ovlap(y,q))
      z->value = intervalHull(z->value,BasicOp::f(p.leftBound(),y));

    y = BasicOp::getRightIntersectionPoint(p);
    if(ovlap(y,q))
      z->value = intervalHull(z->value,BasicOp::f(p.rightBound(),y));
    y = -y; 
    if(ovlap(y,q))
      z->value = intervalHull(z->value,BasicOp::f(p.rightBound(),y));
    //// 
    interval x = BasicOp::getLeftIntersectionPoint(q);
    if(ovlap(x,p))
      z->value = intervalHull(z->value,BasicOp::f(x,q.leftBound()));
    x = -x; 
    if(ovlap(x,p))
      z->value = intervalHull(z->value,BasicOp::f(x,q.leftBound()));

    x = BasicOp::getRightIntersectionPoint(q);
    if(ovlap(x,p))
      z->value = intervalHull(z->value,BasicOp::f(x,q.rightBound()));

    x = -x;
    if(ovlap(x,p))
      z->value = intervalHull(z->value,BasicOp::f(x,q.rightBound()));
    z->isnan = false;
  }  
};

struct AssertEqualityOp{
  inline static void eval(Node* z, Node* l, Node* r) { 
    z->isnan = !intersection(l->value,r->value,z->value); 
  } 
  inline static capd::interval getL(capd::interval l, capd::interval r, capd::interval z) { 
    return getR(l,r,z); 
  }
  inline static capd::interval getR(capd::interval l, capd::interval r, capd::interval z) { 
    if(!intersection(l,r,l)) return r;
    if(!intersection(l,z,l)) return r; 
    return l;
  } 
};

typedef BinaryNode<AddOp> AddNode;
typedef BinaryNode<SubOp> SubNode;
typedef BinaryNode<MulOp> MulNode;
typedef BinaryNode<DivOp> DivNode;
typedef BinaryNode< XByNormOp<XByNorm3Op> > XByNorm3Node;
typedef BinaryNode< XByNormOp<X2ByNorm5Op> > X2ByNorm5Node;
typedef BinaryNode< XByNormOp<XByNorm2Op> > XByNorm2Node;
typedef BinaryNode< NodeWithExplicitBound<SubOp> > SubWithExplicitBoundNode;
typedef BinaryNode<AssertEqualityOp> AssertEqualityNode;

struct SqrNode : public Node{
  using Node::Node;
  void eval(){
    if(!left->isnan){
      this->value = sqr(left->value); 
      this->isnan = false;
    } else
      this->isnan = true;
  }
  Result contract(){
    Result c = NOIMPROVEMENT;
    if(left->value>=0.){
      while(true){
        capd::interval L = sqrt(this->value);
        bool r = L.leftBound()  > left->value.leftBound() or
                 L.rightBound() < left->value.rightBound(); 
        if(!intersection(L,left->value,left->value)) return EMPTY;
        if(r) c=CONTRACTION;
        else break;
      }
    } else if(left->value<=0.){
      while(true){
        capd::interval L = -sqrt(this->value);
        bool r = L.leftBound()  > left->value.leftBound() or
                 L.rightBound() < left->value.rightBound(); 
        if(!intersection(L,left->value,left->value)) return EMPTY;
        if(r) c=CONTRACTION;
        else break;
      }    
    } else {
      while(true){
        capd::interval L = interval(-1,1)*sqrt(this->value);
        bool r = L.leftBound()  > left->value.leftBound() or
                 L.rightBound() < left->value.rightBound(); 
        if(!intersection(L,left->value,left->value)) return EMPTY;
        if(r) c=CONTRACTION;
        else break;
      }    
    }
    Result l = left->contract();
    if(l==EMPTY) return EMPTY;
    return Result(c or l);
  }
};

template<class BasicOp>
struct UnaryNode : public Node{
  using Node::Node;
  void eval(){ 
    if(!left->isnan)
      BasicOp::eval(this,left);
    else 
      this->isnan = true;
  }
  Result contract(){
    Result c = NOIMPROVEMENT;
    while(true){
      capd::interval L = BasicOp::getL(left->value,this->value);
      bool r = L.leftBound()  > left->value.leftBound() or
               L.rightBound() < left->value.rightBound(); 
      if(!intersection(L,left->value,left->value)) return EMPTY;
      if(r) c=CONTRACTION;
      else break;
    }
    Result l = left->contract();
    if(l==EMPTY) return EMPTY;
    return Result(c or l);
  }
};

struct PowOneAndHalfOp{
  inline static void eval(Node* z, Node* l) { 
    capd::interval u = abs(l->value);
    if(!isSingular(u)){
      intersection(power(u,1.5),sqrt(power(u,3)),z->value); 
    } else{ 
      z->value = sqrt(power(u,3)); 
    }
  } 
  inline static capd::interval getL(capd::interval l, capd::interval z) { return power(z,interval(2)/3);} 
};

struct SqrtOp{
  inline static void eval(Node* z, Node* l) { 
    if(!(l->value>=0.)){
      z->isnan = true;
    } else {
      z->value = sqrt(l->value); 
      z->isnan = false;
    }
  }
  inline static capd::interval getL(capd::interval l, capd::interval z) { return sqr(z); } 
};

struct CubeOp{
  inline static void eval(Node* z, Node* l) { z->value = power(l->value,3); }
  inline static capd::interval getL(capd::interval l, capd::interval z) { return power(z,interval(1)/3); } 
};

typedef UnaryNode<PowOneAndHalfOp> PowOneAndHalfNode;
typedef UnaryNode<SqrtOp> SqrtNode;
typedef UnaryNode<CubeOp> CubeNode;

