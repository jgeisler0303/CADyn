// We add a motor torque of 10 Nm on the crank up to 2 s
body[0]->Mext(2)= -4*q[0] - 0.1*qd[0];
