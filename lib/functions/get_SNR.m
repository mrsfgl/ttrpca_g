function [ msqe, psnrval, ssimval, snrval] = get_SNR( X, X0 )

if size(X,ndims(X))==3
    tempX = reshape(X,[],3);
    tempX0 = reshape(X0,[],3);
    msqe = immse(tempX, tempX0);
    for i=1:3
        [psnrval(i), snrval(i)] = psnr(tempX(:,i), tempX0(:,i));
        % trueNorm = norm(X-X0)^2;
        % rel_err = numel(X)*max(abs(X),[],'all')^2 / trueNorm;
        % snr = 10*log10( rel_err );

        % m_x = mean(X);
        % s_x = std(X);
        % m_x0 = mean(X0);
        % s_x0 = std(X0);
        % c = (X-m_x)'*(X0-m_x0);
        % C1 = 1;
        % C2 = 1;
        % ssim = (2*m_x*m_x0)*(2*c+C2)/(m_x^2+m_x0^2+C1)/...
        %     (s_x^2+s_x0^2+C2);
        ssimval(i) = ssim(tempX(:,i), tempX0(:,i));
    end
    psnrval = mean(psnrval);
    snrval = mean(snrval);
    ssimval = mean(ssimval);
else
    msqe = norm(X(:)- X0(:))/norm(X0(:));
    [psnrval, snrval] = psnr(X(:), X0(:));
    ssimval = ssim(X(:), X0(:));
end

end