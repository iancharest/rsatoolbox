% Inverts a colormap.
%
% CW 2015-07
function m_out = invert_colormap(m_in)
    m_out = m_in(end:-1:1, :);
end%function
