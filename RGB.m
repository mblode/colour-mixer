Ta=[5.47813E-05 0.000184722 0.000935514 0.003096265 0.009507714 0.017351596 0.022073595 0.016353161 0.002002407 -0.016177731 -0.033929391 -0.046158952 -0.06381706 -0.083911194 -0.091832385 -0.08258148 -0.052950086 -0.012727224 0.037413037 0.091701812 0.147964686 0.181542886 0.210684154 0.210058081 0.181312094 0.132064724 0.093723787 0.057159281 0.033469657 0.018235464 0.009298756 0.004023687 0.002068643 0.00109484 0.000454231 0.000255925; -4.65552E-05 -0.000157894 -0.000806935 -0.002707449 -0.008477628 -0.016058258 -0.02200529 -0.020027434 -0.011137726 0.003784809 0.022138944 0.038965605 0.063361718 0.095981626 0.126280277 0.148575844 0.149044804 0.14239936 0.122084916 0.09544734 0.067421931 0.035691251 0.01313278 -0.002384996 -0.009409573 -0.009888983 -0.008379513 -0.005606153 -0.003444663 -0.001921041 -0.000995333 -0.000435322 -0.000224537 -0.000118838 -4.93038E-05 -2.77789E-05; 0.00032594 0.001107914 0.005677477 0.01918448 0.060978641 0.121348231 0.184875618 0.208804428 0.197318551 0.147233899 0.091819086 0.046485543 0.022982618 0.00665036 -0.005816014 -0.012450334 -0.015524259 -0.016712927 -0.01570093 -0.013647887 -0.011317812 -0.008077223 -0.005863171 -0.003943485 -0.002490472 -0.001440876 -0.000852895 -0.000458929 -0.000248389 -0.000129773 -6.41985E-05 -2.71982E-05 -1.38913E-05 -7.35203E-06 -3.05024E-06 -1.71858E-06]



for r=0:255
answerstemp=NaN(1,36);
  for g=0:255
    for b=0:255
      saRGB=[r g b]
      spd=ILLSS(Ta,saRGB);
      new=horzcat(spd');
      answerstemp=[answerstemp; new];
    end
   end
dlmwrite(strcat('data/', num2str(sprintf('%3.3d',r))),answerstemp)
end