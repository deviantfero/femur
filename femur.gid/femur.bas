# PROBLEM STAGE 0
#  data: n-nodes n-elements time-delta initial velocity
#  ==============================
*npoin *GenData(End_Time) *GenData(Time_Delta) *GenData(Initial_Velocity) *GenData(Liquid_Density)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# NODES STAGE 1
#  data: node, x, y, z
#  ==============================
*loop nodes
*nodesnum *nodescoord
*end nodes

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# CONECTIVITY STAGE 2
# data: node, c1, c2, c3, c4
# ===============================
*loop elems
*elemsnum *elemsconec(1) *elemsconec(2) *elemsconec(3) *elemsconec(4)
*end elems


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# NOSLIP NODES STAGE 3
# data: nodes
# ===============================
*Set Cond No-Slip *nodes
*loop nodes *OnlyInCond
*nodesnum
*end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# INPUT NODES STAGE 4
# data: node
# ===============================
*Set Cond Input *nodes
*loop nodes *OnlyInCond
*nodesnum
*end

;;;;;;;;;;;;;;;;;
# OUTPUT NODES STAGE 5
# data: node
# ===============================
*Set Cond Output *nodes
*loop nodes *OnlyInCond
*nodesnum
*end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
# NOSLIP ELEMENTS STAGE 6
# data: elements
# ===============================
*Set Cond No-Slip *elems
*loop elems *OnlyInCond
*elemsnum
*end