## Author: Daniele Caratelli
#email: danicaratelli@gmail.com (don't hesitate to contact if you are
#confused, at worst I will not answer you)

#description: Here are the possible transformations of the data

function pc(x);
    return y = 100*(x[2:end]-x[1:end-1])./x[1:end-1];
end    

function pcq(x);
    return y = 400*(x[2:end]-x[1:end-1])./x[1:end-1];
end

function lin(x);
    return y = x;
end

function ln(x);
    return y = log(x);
end
