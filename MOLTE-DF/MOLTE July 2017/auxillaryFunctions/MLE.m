function [ alpha,betas,ave ] = MLE( a,b,f )

    if isempty(b) % one dimension
        [t1,t2,ave]=MLE1d(a,f);
        alpha=t1;
        betas=t2;
    else   % two dimensions
        [t1,t2,t3,ave]=MLE2d(a,b,f);
        alpha=t1;
        betas=[t2 t3];
    end
    
    if alpha<1e-10
        alpha=1;
    end
end

