function v=BSCallVol(p,S,K,v0,r,t)

obj=@(iv)(BSCall(S,K,iv,r,t)-p)^2;

v = fminsearch(obj,v0);