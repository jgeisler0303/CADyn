// We add a motor torque of 10 Nm on the crank up to 2 s
body[0]->Fext(0)= -q[0] -0.01*qd[0];
