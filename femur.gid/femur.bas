# PROBLEM
#  data: n-nodes n-elements time-delta initial velocity
#  ==============================
*npoin *nelem *GenData(Time_Delta) *GenData(Initial_Velocity) *GenData(Liquid_Densitypp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# NODES
#  data: node, x, y, z
#  ==============================
*loop nodes
*nodesnum *nodescoord
*end nodes

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# CONECTIVITY
# data: node, c1, c2, c3, c4
# ===============================
*loop elems
*elemsnum *elemsconec(1) *elemsconec(2) *elemsconec(3) *elemsconec(4)
*end elems

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# NOSLIP NODES
# data: node
# ===============================
*Set Cond No-Slip *nodes
*loop nodes *OnlyInCond
*nodesnum
*end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# INPUT NODES
# data: node
# ===============================
*Set Cond Input *nodes
*loop nodes *OnlyInCond
*nodesnum
*end

;;;;;;;;;;;;;;;;;
# OUTPUT NODES
# data: node
# ===============================
*Set Cond Output *nodes
*loop nodes *OnlyInCond
*nodesnum
*end
