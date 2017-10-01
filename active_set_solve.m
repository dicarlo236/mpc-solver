%some epsilons
epsilon_constraint = 0.00000001;
epsilon_lagrange =   0.00000001;

%load up some example data (including a solution) from cheetah simulation.
%this is in the form:
%min 0.5*x'*H*x + x'*g
%subject to lb < A*x < ub
load('example_data.mat');

%make equality constraints inequality constraints with small bounds
s.ub(s.ub==0)=.00002;

%create active set, start with all constraints inactive (W = 0)
W_upper = zeros(size(s.A,1),1);
W_lower = zeros(size(W_upper,1),1);

%If desired, initialize W_upper and W_lower to something more interesting.

%W will be -1 for a lower constraint, 0 for no constraint, 1 for an upper
%constraint, and -1 for both constraints.
W = -W_lower + W_upper - and(W_lower, W_upper);

cindex = 1:180;

%current solution - starting with something feasible
old_U = repmat([0.00001 0.00001 0.00005]',108/3,1);

tic
for i = 1:2000
    
   %disp(['ITERATION: ' num2str(i)]);
   
   %using the active set of constraints, construct the equality
   %constrainted problem:
   %min 0.5*x'*H*x + g'*x
   %subject to: Ex = d
   E = [s.A(W == 1,:);s.A(W == -1,:)];
   d = [s.ub(W == 1);s.lb(W == -1)];
   %a map between constraints in E and constraints in A.
   EA_index_l = cindex(W==-1);
   EA_index_u = cindex(W==1);
   
   %solve the equality constrained problem
   KKT = [s.H E'; E zeros(size(E,1), size(s.H,1)+size(E,1)-size(E,2)) ];
   GD = [-s.g; d];
   UL = (KKT)\GD;

   %find the direction of the solution
   direction = UL(1:108) - old_U;
   
   %Lets see if this solution is still within constraints:
   lb_check = min(s.A * UL(1:108) - s.lb);
   ub_check = min(s.ub - s.A * UL(1:108));
   if min(lb_check, ub_check) > -epsilon_constraint
       %if so, use it to update solution
       old_U = UL(1:108);
       
       %we can try to remove some constraints
       %disp('CHECK REMOVE CONSTRAINT');
       num_remove = 0;
       %look at the lagrange multiplier to see if moving off the constraint
       %will help us
       lagranges = UL(109:end);
       
       %first check the upper bounds:
       for l = 1:length(EA_index_u)
          if(-lagranges(l) > epsilon_lagrange)
              num_remove = num_remove + 1;
              W(EA_index_u(l)) = 0;
          end
       end
       
       %check the lower bounds:
       for l = 1:length(EA_index_l)
           if(-lagranges(length(EA_index_u) + l) < -epsilon_lagrange)
               num_remove = num_remove + 1;
               W(EA_index_l(l)) = 0;
           end
       end
       
       
       if(num_remove == 0)
           toc
           %If we didn't remove anything, we're all done!
           disp('ACCURACY');
           max(abs([old_U(1:108) - s.result]))
           error('DONE!');
       end
       
   %If the solution to the equality constrained equation pushed us over 
   %a constraint
   else
       %disp('SCALE!');
   %now we can scale the direction such that we hit exactly one new
   %constraint.  First, the lower bounds:
   
   %Only check constraints we're moving toward, startin with lower bounds:
   moving_toward_l = and(s.A * direction < 0, W ~= -1);
   if(sum(moving_toward_l) > 0)
   A_sl = s.A(moving_toward_l,:);
   lb_sl = s.lb(moving_toward_l);
   [scale_required_l,l_index] = min((A_sl * old_U - lb_sl)./(-A_sl*direction));
   MA_index_l = cindex(moving_toward_l);
   end
   
   %upper bounds
   moving_toward_u = and(s.A * direction > 0, W ~= 1);
   if(sum(moving_toward_u) > 0)
   A_sl = s.A(moving_toward_u,:);
   ub_sl = s.ub(moving_toward_u,:);
   [scale_required_u,u_index] = min((ub_sl - A_sl*old_U)./(A_sl*direction));
   MA_index_u = cindex(moving_toward_u);
   end
   
   %pick the lowest scale required to hit one constraint exactly
   scale = min(scale_required_l,scale_required_u);
   
   %disp(['scale: ' num2str(scale)]);
   
   %disp(['sum W: ' num2str(sum(W))]);

   %scale the answer
   old_U = old_U(1:108) + scale .* direction;
   
   %add the constraint
   %if we're driven by a lower bound:
   if scale_required_l < scale_required_u
       %disp('add l');
       W(MA_index_l(l_index)) = -1;
   else
       %disp('add u');
       W(MA_index_u(u_index)) = 1;
   end
   
   end
end