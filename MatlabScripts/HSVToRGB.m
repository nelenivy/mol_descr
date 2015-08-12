function [R, G, B] = HSVToRGB(H, S, V) 
	i = floor(H * 6 / 256); % 0..5
	f = H * 6 - i*256; % 0..255
	p = floor(V * (255 - S)/255); % 0..255
	q = floor(V * (65280 - f * S)/65280); % 65280 = 255*256
	t = floor(V * (65280 - (256 - f) * S)/65280);

    reminder = floor(mod(i, 6));
    
    if (reminder == 0)
      R = V; G = t; B = p;
    elseif (reminder == 1)
         R = q; G = V; B = p;
    elseif (reminder == 2)
        R = p; G = V; B = t;
    elseif (reminder == 3)
        R = p; G = q; B = V;
    elseif (reminder == 4)
        R = t; G = p; B = V;
    else %if (reminder == 5)
        R = V; G = p; B = q;
    end
    
    if (G ~=  43.1373)
        G;
    end
    
     if (B ~=  43.1373)
        B;
    end
    
    R = max(0.0, min(R/ 255.0, 1.0));
    G = max(0.0, min(G/ 255.0, 1.0));
    B = max(0.0, min(B/ 255.0, 1.0));
   end
	