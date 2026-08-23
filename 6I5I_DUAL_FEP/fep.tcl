# fep.tcl -- runFEP helper for NAMD alchemical FEP (dual topology)
# 16 lambda values (15 windows), denser near the endpoints where dE/dlambda
# changes fastest (electrostatics decoupling).
set fep_lambdas {0.0 0.045 0.09 0.14546 0.22425 0.30303 0.38182 0.46061 \
                 0.5394 0.61819 0.697 0.77576 0.85455 0.91 0.955 1.0}

proc runFEP {nsteps {direction forward}} {
    global fep_lambdas
    set n [llength $fep_lambdas]
    for {set i 0} {$i < $n - 1} {incr i} {
        if {$direction eq "backward"} {
            set lambda1 [lindex $fep_lambdas [expr {$n - 1 - $i}]]
            set lambda2 [lindex $fep_lambdas [expr {$n - 2 - $i}]]
        } else {
            set lambda1 [lindex $fep_lambdas $i]
            set lambda2 [lindex $fep_lambdas [expr {$i + 1}]]
        }
        alchLambda  $lambda1
        alchLambda2 $lambda2
        print "FEP window $i: lambda $lambda1 -> $lambda2"
        run $nsteps
    }
}
