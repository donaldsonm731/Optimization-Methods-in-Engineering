%% Problem 5
w = [19 15 20 8 5 7 3 2 4];
v = [380 225 320 96 70 126 30 22 68];
value_star = 0;

for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                for m = 0:1
                    for n = 0:1
                        for o = 0:1
                            for p = 0:1
                                for q = 0:1
                                    weight = sum(w.*[i,j,k,l,m,n,o,p,q]);
                                    value = sum(v.*[i,j,k,l,m,n,o,p,q]);

                                    if weight < 40 && value > value_star
                                        value_star = value;
                                        x_star = [i,j,k,l,m,n,o,p,q];
                                        weight_star = weight;

                                    end
                                     

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
