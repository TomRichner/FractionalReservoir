 I talked to a collaborator who is trying to study the effect of SFA & STD and fractional-order dynamics in reservoir computing.  I'd like to create a echo
  state network which shows that networks with multiple timescale SFA & STD or a fractional-order network (e.g. using the GL fractional derivative) can do
  things that conventional RNNs cannot.  For instance, fractional order derivatives typically have a phase advance relative to sinusoidal inputs, irrespective
  of the input frequency.  In comparison, an LTI system has a fixed set of eigenvalues, each with an associated frequency.  So, if we send in a cyclic
  stimulus with different levels of temporal stretching, an LTI network should respond differently, depending on the time stretch and essentially have
  different encoding for each.  I'd expect a network composed of GL fractional-order neurons would potentially have more temporal invariance and therefore as
  the input is time stretched, the output would also be time stretched, but the encoding would be fairly invariant.  To test this, I'd like to set up an ESN
  in which we train at one time scale, but then stretch the input and test a different time scale.