# i n c l u d e   " c u t s . h " 
 u s i n g   n a m e s p a c e   R o o F i t ; 
 
 v o i d   S i g n a l _ m b c _ f i t ( v o i d ) { 
     T F i l e   * i f i l e   =   T F i l e : : O p e n ( " / h o m e / v i t a l y / B 0 t o D h 0 / T M V A / f i l _ b 2 d h _ s i g . r o o t " ) ; 
     T T r e e   * t r e e   =   ( T T r e e * ) i f i l e - > G e t ( " T E v e n t " ) ; 
 
     R o o C a t e g o r y   b 0 f ( " b 0 f " , " b 0 f " ) ; 
     b 0 f . d e f i n e T y p e ( " g o o d " , 1 ) ; 
     b 0 f . d e f i n e T y p e ( " b a d " , 5 ) ; 
     b 0 f . d e f i n e T y p e ( " f s r " , 1 0 ) ; 
 
     R o o A r g S e t   a r g s e t ; 
     R o o R e a l V a r   m b c ( " m b c " , " M _ { b c } " , 5 . 2 0 , 5 . 3 0 , " G e V " ) ;   a r g s e t . a d d ( m b c ) ; 
     R o o R e a l V a r   d e ( " d e " , " # D e l t a E " , - 0 . 3 , 0 . 3 , " G e V " ) ; 
     a r g s e t . a d d ( d e ) ; 
 / /       d e . s e t R a n g e ( " s i g n a l " , d e _ m i n , d e _ m a x ) ; 
 / /       d e . s e t R a n g e ( " f i t " , - 0 . 1 5 , 0 . 3 ) ; 
     R o o R e a l V a r   m d ( " m d " , " m d " , D M a s s - m d _ c u t , D M a s s + m d _ c u t , " G e V " ) ;   a r g s e t . a d d ( m d ) ; 
     R o o R e a l V a r   m k ( " m k " , " m k " , K M a s s - m k _ c u t , K M a s s + m k _ c u t , " G e V " ) ;   a r g s e t . a d d ( m k ) ; 
     R o o R e a l V a r   m p i 0 ( " m p i 0 " , " m p i 0 " , P i 0 M a s s - m p i 0 _ c u t , P i 0 M a s s + m p i 0 _ c u t , " G e V " ) ;   a r g s e t . a d d ( m p i 0 ) ; 
     R o o R e a l V a r   b d t g s ( " b d t g s " , " b d t g s " , 0 . 9 8 , 1 . ) ;   a r g s e t . a d d ( b d t g s ) ; 
     R o o R e a l V a r   a t c k p i _ m a x ( " a t c k p i _ m a x " , " a t c k p i _ m a x " , 0 . , a t c k p i _ c u t ) ;   a r g s e t . a d d ( a t c k p i _ m a x ) ; 
 
     a r g s e t . a d d ( b 0 f ) ; 
 
     R o o D a t a S e t   d s ( " d s " , " d s " , t r e e , a r g s e t ) ; 
     d s . P r i n t ( ) ; 
 
 / /     g R O O T - > P r o c e s s L i n e ( " . L   p d f s / R o o S t u d e n t s G a u s s 1 D . c x x + " ) ; 
 
     R o o R e a l V a r   m b c 0 ( " m b c 0 " , " m b c 0 " , 5 . 2 8 6 , 5 . 2 6 , 5 . 3 0 ) ; 
 / /       R o o R e a l V a r   w ( " w " , " w " , 0 . 0 2 1 , 0 . , 0 . 5 ) ; 
 / /       R o o R e a l V a r   d w ( " d w " , " d w " , 0 . 0 , - 0 . 0 1 , 0 . 0 1 ) ;   d w . s e t C o n s t a n t ( k T R U E ) ; 
 / /       R o o R e a l V a r   n h ( " n h " , " n h " , 2 . , 0 . 1 , 1 0 . ) ; 
 / /       R o o R e a l V a r   n l ( " n l " , " n l " , 2 . , 0 . 1 , 1 0 . ) ; 
 / /       R o o R e a l V a r   f ( " f " , " f " , 0 . 3 , 0 . , 1 . ) ; 
 / /       R o o R e a l V a r   d m b c 0 ( " d m b c 0 " , " d m b c 0 " , 0 . , - 0 . 1 , 0 . 1 ) ;   d m b c 0 . s e t C o n s t a n t ( k T R U E ) ; 
 / /       R o o R e a l V a r   s g ( " s g " , " s g " , 1 . , 0 . 9 , 1 . 1 ) ;   s g . s e t C o n s t a n t ( k T R U E ) ; 
     R o o R e a l V a r   s l ( " s l " , " s l " , 0 . 0 1 3 5 , 0 . , 0 . 5 ) ; 
     R o o R e a l V a r   s r ( " s r " , " s r " , 0 . 0 0 1 , 0 . , 0 . 5 ) ; 
     R o o B i f u r G a u s s   b g ( " b g " , " b g " , m b c , m b c 0 , s l , s r ) ; 
 
 / /     R o o S t u d e n t s G a u s s 1 D   p d f ( " p d f " , " p d f " , m b c , m b c 0 , w , d w , n h , n l , f , d m b c 0 , s g ) ; 
 
     R o o R e a l V a r   m b c 0 0 ( " m b c 0 0 " , " m b c 0 0 " , 5 . 2 8 0 , 5 . 2 6 , 5 . 3 0 ) ; 
     R o o R e a l V a r   s l l ( " s l l " , " s l l " , 0 . 0 0 8 7 , 0 . , 0 . 5 ) ; 
     R o o R e a l V a r   s r r ( " s r r " , " s r r " , 0 . 0 0 2 8 , 0 . , 0 . 5 ) ; 
     R o o B i f u r G a u s s   b g g ( " b g g " , " b g g " , m b c , m b c 0 0 , s l l , s r r ) ; 
   
 / /       R o o R e a l V a r   m b c C B l ( " m b c C B l " , " m b c C B l " , 5 . 2 8 , 5 . 2 7 , 5 . 2 9 ) ; 
 / /       R o o R e a l V a r   s C B l ( " s C B l " , " s C B l " , 0 . 0 3 2 , 0 . , 0 . 5 ) ; 
 / /       R o o R e a l V a r   n l ( " n l " , " n l " , 1 2 . 1 , 0 . , 1 0 0 . ) ; 
 / /       R o o R e a l V a r   a l p h a l ( " a l p h a l " , " a l p h a l " , 0 . 7 2 , - 1 0 . , 1 0 . ) ; 
 / /   
 / /       R o o R e a l V a r   m b c C B r ( " m b c C B r " , " m b c C B r " , 5 . 2 8 , 5 . 2 7 , 5 . 2 9 ) ; 
 / /       R o o R e a l V a r   s C B r ( " s C B r " , " s C B r " , 0 . 0 6 6 , 0 . , 0 . 5 ) ; 
 / /       R o o R e a l V a r   n r ( " n r " , " n r " , 1 9 . 7 , 0 . , 1 0 0 . ) ; 
 / /       R o o R e a l V a r   a l p h a r ( " a l p h a r " , " a l p h a r " , - 0 . 7 5 , - 1 0 . , 1 0 . ) ; 
 / /       
 / /       R o o C B S h a p e   C B l ( " C B l " , " C B l " , m b c , m b c C B l , s C B l , a l p h a l , n l ) ; 
 / /       R o o C B S h a p e   C B r ( " C B r " , " C B r " , m b c , m b c C B r , s C B r , a l p h a r , n r ) ; 
 / /       
 / /       R o o R e a l V a r   f C B l ( " f C B l " , " f C B l " , 0 . 2 0 , 0 . , 1 . ) ; 
 / /       R o o R e a l V a r   f C B r ( " f C B r " , " f C B r " , 0 . 2 0 , 0 . , 1 . ) ; 
 
       R o o R e a l V a r   f ( " f " , " f " , 0 . 4 9 7 , 0 . , 1 . ) ;       
       R o o A d d P d f   p d f ( " p d f " , " p d f " , R o o A r g L i s t ( b g g , b g ) , R o o A r g S e t ( f ) ) ; 
 
 / /       R o o R e a l V a r   m e a n ( " m e a n " , " m e a n " , 5 . 2 8 , 5 . 2 7 , 5 . 2 9 ) ; 
 / /       R o o R e a l V a r   s i g m a ( " s i g m a " , " s i g m a " , 0 . 0 0 5 , 0 . 0 0 1 , 0 . 5 ) ; 
 / /       R o o G a u s s i a n   p d f ( " p d f " , " p d f " , m b c , m e a n , s i g m a ) ; 
     
       s t r i n g s t r e a m   o u t ;     
 / /       o u t . s t r ( " " ) ; 
 / /       o u t   < <   " d e < "   < <   d e _ m a x   < <   " & & d e > "   < <   d e _ m i n ; 
 / /       c o n s t   i n t   g e n _ b a c k   =   d s . s u m E n t r i e s ( o u t . s t r ( ) . c _ s t r ( ) ) ; 
 / /       c o u t   < <   " G e n   b a c k g r o u n d :   "   < <   g e n _ b a c k   < <   e n d l ; 
     
 / /     p d f . f i t T o ( d s , R a n g e ( " f i t " ) , V e r b o s e ( ) , T i m e r ( t r u e ) ) ; 
     p d f . f i t T o ( d s , V e r b o s e ( ) , T i m e r ( t r u e ) ) ; 
     
     / / / / / / / / / / / / / 
     / /     P l o t s     / / 
     / / / / / / / / / / / / / 
     R o o P l o t *   d e F r a m e   =   m b c . f r a m e ( ) ; 
     d s . p l o t O n ( d e F r a m e , D a t a E r r o r ( R o o A b s D a t a : : S u m W 2 ) , M a r k e r S i z e ( 1 ) ) ; 
     p d f . p l o t O n ( d e F r a m e , C o m p o n e n t s ( b g g ) , L i n e W i d t h ( 1 ) , L i n e S t y l e ( k D a s h e d ) ) ; 
     p d f . p l o t O n ( d e F r a m e , C o m p o n e n t s ( b g ) , L i n e W i d t h ( 1 ) , L i n e S t y l e ( k D a s h e d ) ) ; 
     p d f . p l o t O n ( d e F r a m e , L i n e W i d t h ( 2 ) ) ; 
 
     R o o H i s t *   h d e p u l l   =   d e F r a m e - > p u l l H i s t ( ) ; 
     R o o P l o t *   d e P u l l   =   m b c . f r a m e ( T i t l e ( " # D e l t a   E   p u l l   d i s t r i b u t i o n " ) ) ; 
     d e P u l l - > a d d P l o t a b l e ( h d e p u l l , " P " ) ; 
     d e P u l l - > G e t Y a x i s ( ) - > S e t R a n g e U s e r ( - 5 , 5 ) ; 
 
     T C a n v a s *   c m   =   n e w   T C a n v a s ( " # D e l t a   E ,   S i g n a l " , " # D e l t a   E ,   S i g n a l " , 6 0 0 , 7 0 0 ) ; 
     c m - > c d ( ) ; 
 
     T P a d   * p a d 3   =   n e w   T P a d ( " p a d 3 " , " p a d 3 " , 0 . 0 1 , 0 . 2 0 , 0 . 9 9 , 0 . 9 9 ) ; 
     T P a d   * p a d 4   =   n e w   T P a d ( " p a d 4 " , " p a d 4 " , 0 . 0 1 , 0 . 0 1 , 0 . 9 9 , 0 . 2 0 ) ; 
     p a d 3 - > D r a w ( ) ; 
     p a d 4 - > D r a w ( ) ; 
 
     p a d 3 - > c d ( ) ; 
     p a d 3 - > S e t L e f t M a r g i n ( 0 . 1 5 ) ; 
     p a d 3 - > S e t F i l l C o l o r ( 0 ) ; 
 
     d e F r a m e - > G e t X a x i s ( ) - > S e t T i t l e S i z e ( 0 . 0 5 ) ; 
     d e F r a m e - > G e t X a x i s ( ) - > S e t T i t l e O f f s e t ( 0 . 8 5 ) ; 
     d e F r a m e - > G e t X a x i s ( ) - > S e t L a b e l S i z e ( 0 . 0 4 ) ; 
     d e F r a m e - > G e t Y a x i s ( ) - > S e t T i t l e O f f s e t ( 1 . 6 ) ; 
     d e F r a m e - > D r a w ( ) ; 
 
     T P a v e T e x t   * p t   =   n e w   T P a v e T e x t ( 0 . 5 5 , 0 . 8 8 , 0 . 9 5 , 0 . 9 4 , " b r N D C " ) ; 
     p t - > S e t F i l l C o l o r ( 0 ) ; 
     p t - > S e t T e x t A l i g n ( 1 2 ) ; 
     o u t . s t r ( " " ) ; 
     o u t   < <   " # c h i ^ { 2 } / n . d . f   =   "   < <   d e F r a m e - > c h i S q u a r e ( ) ; 
     p t - > A d d T e x t ( o u t . s t r ( ) . c _ s t r ( ) ) ; 
 / /       o u t . s t r ( " " ) ; 
 / /       o u t   < <   " S i g n a l   r e g i o n :   ( "   < <   d e _ m i n   < <   " , "   < <   d e _ m a x   < <   " ) " ; 
 / /       p t - > A d d T e x t ( o u t . s t r ( ) . c _ s t r ( ) ) ; 
 / /       o u t . s t r ( " " ) ; 
 / /       o u t   < <   " E v e n t s   i n   t h e   S R :   "   < <   g e n _ b a c k ; 
 / /       p t - > A d d T e x t ( o u t . s t r ( ) . c _ s t r ( ) ) ; 
     p t - > D r a w ( ) ; 
 
 / /       T L i n e   * d e _ l i n e L E F T   =   n e w   T L i n e ( d e _ m i n , 0 , d e _ m i n , 0 . 3 * g e n _ b a c k ) ; 
 / /   / /     T L i n e   * d e _ l i n e L E F T   =   n e w   T L i n e ( d e _ m i n , 0 , d e _ m i n , 1 8 ) ; 
 / /       d e _ l i n e L E F T - > S e t L i n e C o l o r ( k R e d ) ; 
 / /       d e _ l i n e L E F T - > S e t L i n e S t y l e ( 1 ) ; 
 / /       d e _ l i n e L E F T - > D r a w ( ) ; 
 / /       
 / /       T L i n e   * d e _ l i n e R I G H T   =   n e w   T L i n e ( d e _ m a x , 0 , d e _ m a x , 0 . 3 * g e n _ b a c k ) ; 
 / /   / /     T L i n e   * d e _ l i n e R I G H T   =   n e w   T L i n e ( d e _ m a x , 0 , d e _ m a x , 1 8 ) ; 
 / /       d e _ l i n e R I G H T - > S e t L i n e C o l o r ( k R e d ) ; 
 / /       d e _ l i n e R I G H T - > S e t L i n e S t y l e ( 1 ) ; 
 / /       d e _ l i n e R I G H T - > D r a w ( ) ; 
     
     p a d 4 - > c d ( ) ;   p a d 4 - > S e t L e f t M a r g i n ( 0 . 1 5 ) ;   p a d 4 - > S e t F i l l C o l o r ( 0 ) ; 
     d e P u l l - > S e t M a r k e r S i z e ( 0 . 0 5 ) ;   d e P u l l - > D r a w ( ) ; 
     T L i n e   * d e _ l i n e U P   =   n e w   T L i n e ( 5 . 2 , 3 , 5 . 2 9 , 3 ) ; 
     d e _ l i n e U P - > S e t L i n e C o l o r ( k B l u e ) ; 
     d e _ l i n e U P - > S e t L i n e S t y l e ( 2 ) ; 
     d e _ l i n e U P - > D r a w ( ) ; 
     T L i n e   * d e _ l i n e   =   n e w   T L i n e ( 5 . 2 , 0 , 5 . 2 9 , 0 ) ; 
     d e _ l i n e - > S e t L i n e C o l o r ( k B l u e ) ; 
     d e _ l i n e - > S e t L i n e S t y l e ( 1 ) ; 
     d e _ l i n e - > S e t L i n e W i d t h ( ( W i d t h _ t ) 2 . ) ; 
     d e _ l i n e - > D r a w ( ) ; 
     T L i n e   * d e _ l i n e D O W N   =   n e w   T L i n e ( 5 . 2 , - 3 , 5 . 2 9 	 , - 3 ) ; 
     d e _ l i n e D O W N - > S e t L i n e C o l o r ( k B l u e ) ; 
     d e _ l i n e D O W N - > S e t L i n e S t y l e ( 2 ) ; 
     d e _ l i n e D O W N - > D r a w ( ) ; 
 
     c m - > U p d a t e ( ) ; 
 } 
 