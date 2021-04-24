pick=10:60:300
%% 作阶数-误差关系图
if(1)
    figure(3);
    for i=1:length(pick)
        S3=MinOrder:Internal:MaxOrder;
        %S3(Err==Minerr)
        plot(S3,10*log10(Err(pick(i),:)),'+-');
        xlabel('order'); ylabel('rmserr(dB)'); 
        grid on;
        hold on;
    end
end