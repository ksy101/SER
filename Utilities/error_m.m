function em = error_m(image_a, image_b)
em = norm(image_a - image_b)/norm(image_b);