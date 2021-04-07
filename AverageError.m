function MAE = AverageError (time, expdata, modeldata, it)
    n=length(time);
    time = time*it;
    MAE = 0;
    for i = 1:n
        MAE = MAE+ abs (expdata(i) - modeldata(time(i)+1));
    end
    MAE = MAE/n;
end