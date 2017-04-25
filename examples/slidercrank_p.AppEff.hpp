    // We add a motor torque of 10 Nm on the crank up to 2 s
    if (t<2.0) body[0]->Mext(2)=10.0; else body[0]->Mext(2)=0.0;
