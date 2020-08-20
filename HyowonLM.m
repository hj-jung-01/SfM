function [x_opt, Residual, J_LM] = HyowonLM(func, x, params)

%     HyowonLM : non-linear least squares problem solver
%     E-mail : hhw8803@gmail.com

%     Copyright (C) 2015, Hyowon Ha
%     All rights reserved. 
%  
%     Redistribution and use in source and binary forms, with or without  
%     modification, are permitted provided that the following conditions are  
%     met: 
%  
%         * Redistributions of source code must retain the above copyright  
%         notice, this list of conditions and the following disclaimer. 
%         * Redistributions in binary form must reproduce the above copyright  
%         notice, this list of conditions and the following disclaimer in the  
%         documentation and/or other materials provided with the distribution. 
%         * Neither the name of the <ORGANIZATION> nor the names of its  
%         contributors may be used to endorse or promote products derived from  
%         this software without specific prior written permission. 
%  
%     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS  
%     IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED  
%     TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A  
%     PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER  
%     OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  
%     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,  
%     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  
%     PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
%     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  
%     NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  
%     SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

l2_p=norm(x);

iter = 0;
max_iter = 100;
lambda = 0.01;
tolX = 1e-6;
tolFun = 1e-6;
tolOpt = 1e-4 * tolFun;
sqrtEps=sqrt(eps);

[F,J] = feval(@(t) func(t, params), x); % target function
sqSumF = dot(F,F);
disp([num2str(iter) ' : ' num2str(sqSumF)]);
Residual = sqSumF;
H = J'*J;
JtF = J'*F;
linf_JtF=max(norm(JtF,Inf),sqrtEps); 

while iter < max_iter
    iter = iter+1;
    
    H_LM = H + sparse(1:length(x),1:length(x),lambda,length(x),length(x));
    dp = -H_LM \ JtF;
    x_LM = x+dp;
    l2_dp=norm(dp);
    
    if ( l2_dp < tolX*(sqrtEps + l2_p) )
        disp('Finished (tolX)');
        break;
    end
    
    [F_LM,J_LM] = feval(@(t) func(t, params), x_LM);
    sqSumF_LM = dot(F_LM,F_LM);
    JtF_LM = J_LM'*F_LM;
    linf_JtF_LM=max(norm(JtF_LM,Inf),sqrtEps);
    
    if (linf_JtF_LM < tolOpt * linf_JtF)
        disp('Finished (tolOpt)');
        break;
    end
    
    if ( abs(sqSumF_LM - sqSumF) <= tolFun*sqSumF )
        disp('Finished (tolFun)');
        break;
    end
    
    if (sqSumF_LM < sqSumF)
        lambda=lambda*0.1;
        x=x_LM;
        sqSumF=sqSumF_LM;
        disp([num2str(iter) ' : ' num2str(sqSumF)]);
        Residual = [Residual, sqSumF];
        H = J_LM'*J_LM;
        JtF = JtF_LM;
    else
        lambda=lambda*10;
    end
end

if iter==max_iter
    disp('Finished (max_iter)');
end

x_opt=x_LM;
