function [features] = histFeas(rawdata)

features = [];
for i = 1:4
    data = rawdata(:,i);

    skew = skewness2(data);

    kurt = kurtosis2(data);

    m1 = mean(data);
    m2 = median(data);

    q1 = quantile(data,0.25);
    q2 = quantile(data,0.75);
    iqr = q2-q1;

    sd = std(data);
    variance = var(data);

    energy = sum(data.^2)/length(data);

    temp = [skew, kurt, m1, m2, q1, q2, iqr, sd, variance, energy];
    features = [features, temp];
end

