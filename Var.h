// Wrapper for Node* to simplify user interface
struct Var{
  Var() : node(nullptr){}
  explicit Var(Node* node) : node(node) {}
  explicit Var(capd::interval x) : node(new ConstNode(x)) {
    Node::dag->push_back(node);
  }


  inline capd::interval& x() { return node->value; }
  inline operator interval&() { return node->value;}
  inline void setValue(capd::interval p) { node->setValue(p); }
  Var(const Var&) = default;
  Var(Var&) = default;
  Var(Var&&) = default;
  Var& operator=(const Var&) = default;
  Var& operator=(Var&&) = default;

  Node* node;
};

inline bool operator<(const Var& x, const Var& y){
  return x.node->value < y.node->value;
}

inline bool operator>(const Var& x, const Var& y){
  return x.node->value > y.node->value;
}

template<class NodeT>
Var binary_op(Var left, Var right){
  Var r = Var(new NodeT(left.node,right.node));
  Node::dag->push_back(r.node);
  return r;  
}

template<class NodeT>
Var unary_op(Var left){
  Var r = Var(new NodeT(left.node,nullptr));
  Node::dag->push_back(r.node);
  return r;  
}

Var operator+(Var left, Var right){
  return binary_op<AddNode>(left,right);
}
Var& operator+=(Var& left, Var right){
  left = left + right;
  return left;
}

Var operator-(Var left, Var right){
  return binary_op<SubNode>(left,right);
}

Var& operator-=(Var& left, Var right){
  left = left - right;
  return left;
}
Var operator*(Var left, Var right){
  return binary_op<MulNode>(left,right);
}
Var operator/(Var left, Var right){
  return binary_op<DivNode>(left,right);
}
Var sqr(Var left){
  return unary_op<SqrNode>(left);
}
Var sqrt(Var left){
  return unary_op<SqrtNode>(left);
}
Var cube(Var left){
  return unary_op<CubeNode>(left);
}
Var operator*(capd::interval p, Var right){
  Var left = Var(p);
  return left*right; 
}
Var powOneAndHalf(Var left){
  return unary_op<PowOneAndHalfNode>(left);
}
Var xByNorm3(Var left, Var right){
  return binary_op<XByNorm3Node>(left,right);
}
Var x2ByNorm5(Var left, Var right){
  return binary_op<X2ByNorm5Node>(left,right);
}
Var xByNorm2(Var left, Var right){
  return binary_op<XByNorm2Node>(left,right);
}
Var assertEquality(Var left, Var right){
  return binary_op<AssertEqualityNode>(left,right);
}

Var subWithExplicitBound(Var left, Var right, capd::interval bound){
  SubWithExplicitBoundNode* node = new SubWithExplicitBoundNode(left.node,right.node);
  node->bound = bound;
  Node::dag->push_back(node);
  return Var(node);  
}
