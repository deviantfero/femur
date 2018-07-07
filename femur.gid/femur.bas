#  n-points   n-nodes  n-elements
#  ======================================
   *npoin      *nnode  *nelem

# ========================================
# n      | x      | y        |    z
# ========================================
*loop nodes
  *nodesnum *nodescoord
*end nodes

# ======= conectivity table ======
*loop elems
  *elemsnum *elemsconec(1) *elemsconec(2) *elemsconec(3) *elemsconec(4)
*end elems
