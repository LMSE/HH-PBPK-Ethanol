function MAE = AveragePercentError (time, expdata, modeldata, it)
    n=length(time);
    time = time*it;
    MAE = 0;
    for i = 1:n
        if (expdata(i) ~= 0 )
        MAE = MAE+ abs (expdata(i) - modeldata(time(i)+1))/expdata(i);
        else 
        MAE = MAE+ abs (expdata(i) - modeldata(time(i)+1))/1;
    end
    MAE = MAE/n;
end